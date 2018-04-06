function [seedPoints, Info] = computeNucleiCenters(BW,options,varargin)
% COMPUTENUCLEICENTERS Find the nuclei centers in a binary image of nuclei.
%
% [seedPoints, Info] = computeNucleiCenters(BW, options)
% [seedPoints, Info] = computeNucleiCenters(BW, options, 'Minimum_Hole_Size', mhs, 'Object_Of_Interest', OoI)
%
% For each object in the image, compute the centers of curvature to use as
% initial points, then find the seed-points for each object.
%
% Input parameters:
% BW : Binary mask of nuclei image.
% options : An instance of class seedPointOptions.
%
% Optional parameter/value pairs:
% 'Minimum_Hole_Size' : The minimum hole size (area) allowed in the mask.
%   Any smaller hole will be filled. Default value is 15.
% 'Object_Of_Interest' : The index of an object of interest. All other
%   objects will be ignored.
%
% Output parameters:
% seedPoints : Nx3 array where the first two columns give the x,y location
%   of a nuclei center and the third column gives the object number to
%   which the nuclei belongs.
% Info : Cell array with a length equal to the number of objects. Each
%   element is the Info structure returned by computeObjectSeedPoints().

% James Kapaldo

[Minimum_Hole_Size, Object_Of_Interest] = parse_inputs(BW,options,varargin{:});

% Get number of rows in image
nRows = size(BW,1);

% Remove small holes from mask
BW = ~bwareaopen(~BW, Minimum_Hole_Size,4);
CC = bwconncomp(BW); % Get connected commponents
pixelList = CC.PixelIdxList;

% Get object scales
objectScale = compute_objectScale(BW, pixelList);

% Compute boundary information
[B,~,K,r0set] = computeBoundaryInformation(BW, objectScale, options);

% If there is an object of interest, remove all the others.
[pixelList, B, K, r0set, objectScale] = getObjectOfInterest(Object_Of_Interest, pixelList, B, K, r0set, objectScale);

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
Que = progres.start();

if Use_Parallel
    % If we are computing in parallel, then first convert the options class
    % element to a structure to prevent reinitiallization on transfer to
    % each worker.
    warning('off','MATLAB:structOnObject')
    options = struct(options);
    warning('on','MATLAB:structOnObject')

    parfor obj = 1:N
        [seedPoints{obj}, Info{obj}] = processObject(pixelList{obj}, nRows, r0set{obj}, useCentroid(obj), obj, options);
        if ~isempty(Que), send(Que, obj), end
    end
else
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

function [Minimum_Hole_Size, Object_Of_Interest] = parse_inputs(BW,options,varargin)

p = inputParser;
p.FunctionName = 'computeNucleiCenters';

    function validate_OoI(t)
        if ~isempty(t) && (numel(t)~=1 || (round(t)-t)>1e-3 || t<=0 || ~isreal(t) || ~isfinite(t))
            error('computeNucleiCenters:badObjectOfInterest','The Object_Of_Interest should be a scalar, real, finite, integer greater than 0.')
        end
    end

addRequired(p,'BW', @(t) validateattributes(t,{'logical'},{'2d'}));
addRequired(p,'options', @(t) isa(t, 'seedPointOptions'))
addParameter(p,'Minimum_Hole_Size', 15, @(t) validateattributes(t,{'double'},{'scalar','nonnegative','real','finite'}))
addParameter(p,'Object_Of_Interest', [], @(t) validate_OoI(t))

parse(p,BW,options,varargin{:})

Minimum_Hole_Size = p.Results.Minimum_Hole_Size;
Object_Of_Interest = p.Results.Object_Of_Interest;

end


%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
