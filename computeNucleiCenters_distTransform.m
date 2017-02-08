function [seedPoints, Info] = computeNucleiCenters_distTransform(I,BW,options)
% COMPUTENUCLEICENTERS Declump the nuclei in an image
%
% [BW,cuts,Info] = computeNucleiCenters(I,BW,options)
%
% I - input image
% BW - object mask for image
% options - declumpOptions class object
%
% BW - object mask after partitioning clumps
% cuts - Nx4 array. cuts(i,1:2) and cuts(i,3:4) give the two vertices of
%        the i'th cut
% Info - cell array of structures giving information

% James Kapaldo



% Get number of rows in image
nRows = size(BW,1);

% Remove small holes from mask
BW = ~bwareaopen(~BW, options.Minimum_Hole_Size,4);
CC = bwconncomp(BW); % Get connected commponents
pixelList = CC.PixelIdxList;

% Get object scales
objectScale = compute_objectScale(BW, pixelList);

% Create confining potential modifier -----------------------------------
% I = imfilter(I, fspecial('gaussian',7,1)); % Smooth the image
% H = homogenize_image(I, BW, CC, options); % Get homogenized image
% H = 1 ./ H; % confining potential
% H = cellfun(@(x) H(x), pixelList,'UniformOutput',false); % Slice potential
H = cell(numel(pixelList),1);

% Compute boundary information
[B,~,K,r0set] = computeBoundaryInformation(BW, objectScale, options);

% If there is an object of interest, remove all the others.
[pixelList, B, K, r0set, H] = getObjectOfInterest(options, pixelList, B, K, r0set, H);

% Use the object centroid as the seed point for any object that is convex
% or smaller than the particle area
area = cellfun(@numel,pixelList); % object areas
isConvex = cellfun(@(x,y) ~any(x > (0.25./y)),K, num2cell(objectScale)); % object convex?
useCentroid = isConvex | (area < pi*options.Wigner_Seitz_Radius.^2);

% Offset the r0set points to coorespond to object origin
objOffset = cellfun(@(x) min(x,[],1,'omitnan') - 1, B,'UniformOutput',false);
r0set = cellfun(@(x,y) x - y, r0set, objOffset,'UniformOutput',false);

% Compute the seed points.
% Note, if you want to also use the seed points to segment the objects, then it would be a good choice to insert the segmentation code in this function "processObjects" and add in any additional input/outputs you need.
[seedPoints, Info] = processObjects(pixelList, H, r0set, useCentroid, nRows, options);

% Offset seed points to image coordinates
seedPoints = cellfun(@(x,y) [x(:,1:2) + y, x(:,3)], seedPoints, objOffset,'UniformOutput',false);

% Offset r0 and r_final to image coordinates
if options.Debug
    for i = 1:numel(Info)
        Info{i}.objOffset = objOffset{i};
    end
end

% Convert seedPoints from cell array to array, Nx3
seedPoints = cat(1,seedPoints{:});

if options.Use_GPU
    gpuDevice([]);
end

if all(cellfun(@isempty, Info))
    Info = [];
end

end
