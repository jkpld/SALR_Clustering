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
nRows = size(I,1);

% Smooth the image
I = imfilter(I, fspecial('gaussian',7,1));

% Remove small holes from mask
BW = ~bwareaopen(~BW, options.Minimum_Hole_Size,4);
CC = bwconncomp(BW); % Get connected commponents
pixelList = CC.PixelIdxList;

% Create interior confining potential -----------------------------------
V_interior = create_base_interior_confining_potential(BW); % Base potential
H = homogenize_image(I, BW, CC, options); % Get homogenized image
V_interior = V_interior ./ H; % confining potential
V_interior = cellfun(@(x) V_interior(x), pixelList,'UniformOutput',false); % Slice potential

% Compute boundary information
[B,~,K,r0set] = computeBoundaryInformation(BW, options);

% If there is an object of interest, remove all the others.
[pixelList, B, K, r0set, V_interior] = getObjectOfInterest(options, pixelList, B, K, r0set, V_interior);

% Use the object centroid as the seed point for any object that is convex
% or smaller than the particle area
area = cellfun(@numel,pixelList); % object areas
isConvex = cellfun(@(x) ~any(x > (0.5/options.Curvature_Max_Radius)),K); % object convex?
useCentroid = isConvex | (area < pi*options.Wigner_Seitz_Radius.^2);

% Offset the r0set points to coorespond to object origin
objOffset = cellfun(@(x) min(x,[],1,'omitnan') - 1, B,'UniformOutput',false);
r0set = cellfun(@(x,y) x - y, r0set, objOffset,'UniformOutput',false);

% Compute the seed points.
% Note, if you want to also use the seed points to segment the objects, then it would be a good choice to insert the segmentation code in this function "processObjects" and add in any additional input/outputs you need.
[seedPoints, Info] = processObjects(pixelList, V_interior, r0set, useCentroid, nRows, options);

% Offset seed points to image coordinates
seedPoints = cellfun(@(x,y) [x(:,1:2) + y, x(:,3)], seedPoints, objOffset,'UniformOutput',false);

% Offset r0 and r_final to image coordinates
if options.Debug
    for i = 1:numel(Info)
        Info{i}.r0 = Info{i}.r0 + objOffset{i};
        Info{i}.r_final = Info{i}.r_final + objOffset{i};
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
