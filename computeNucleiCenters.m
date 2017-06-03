function [seedPoints, Info] = computeNucleiCenters(I,BW,options)
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

% Compute boundary information
[B,~,K,r0set] = computeBoundaryInformation(BW, objectScale, options);

% If there is an object of interest, remove all the others.
[pixelList, B, K, r0set, objectScale] = getObjectOfInterest(options, pixelList, B, K, r0set, objectScale);

% Use the object centroid as the seed point for any object that is convex
% or smaller than the particle area
area = cellfun(@numel,pixelList); % object areas
isConvex = cellfun(@(x,y) ~any(x(~isnan(x)) > (0.25./y)),K, num2cell(objectScale)); % object convex?
useCentroid = isConvex | (area < pi*options.Wigner_Seitz_Radius.^2);

% Offset the r0set points to coorespond to object origin
objOffset = cellfun(@(x) min(x,[],1,'omitnan') - 1, B,'UniformOutput',false);
r0set = cellfun(@(x,y) x - y, r0set, objOffset,'UniformOutput',false);

% Compute the seed points.
N = numel(pixelList);
verbose = options.Verbose;
Use_Parallel = options.Use_Parallel;
if N == 1
    Use_Parallel = 0;
end

% Initialize sliced variables
Info = cell(N,1);
seedPoints = cell(N,1);

% Create a progress monitor
progres = displayProgress(N, 10, verbose, Use_Parallel);

if Use_Parallel
    % If we are computing in parallel, then first convert the options class
    % element to a structure to prevent reinitiallization on transfer to
    % each worker.
    warning('off','MATLAB:structOnObject')
    options = struct(options);
    warning('on','MATLAB:structOnObject')

    Que = progres.start();    
    parfor obj = 1:N
        [seedPoints{obj}, Info{obj}] = processObject(pixelList{obj}, nRows, r0set{obj}, useCentroid(obj), obj, options);
        if ~isempty(Que), send(Que, obj), end
    end    
else
    progres.start();
    for obj = 1:N
        [seedPoints{obj}, Info{obj}] = processObject(pixelList{obj}, nRows, r0set{obj}, useCentroid(obj), obj, options);
        progres.iteration_end()
    end
end

delete(progres)

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

function [seedPoints, Info] = processObject(pixList, nRows, r0set, useCentroid, obj, options)

    % Create mask and interior potential images for the object
    objBW = createObjectImages(pixList, nRows, true(numel(pixList),1));
    [seedPoints, Info] = computeObjectSeedPoints(objBW, options, 'r0set', r0set, 'useCentroid', useCentroid, 'objNumber', obj);

    % Add object number as third column. (This is mostly just helpful when comparing against truth data, as the truth data is labeled by each object.)
    seedPoints = [seedPoints, obj*ones(size(seedPoints,1),1)];

    % ====================================================================
    % This could be a good place to put object segmentation code. Just add
    % any additional inputs/outputs you need.
    % ====================================================================
end
