function [dV, r0, problem_scales, V, Info] = setup_problem(binned_data, data_limits, options, r0set)
% SETUP_PROBLEM Generate the confining force, the initial particle
% locations, and the scaling factors for the simulation.
%
% [dV, r0, problem_scales, V, Info] = setup_problem(binned_data, data_range, options)
% [dV, r0, problem_scales, V, Info] = setup_problem(binned_data, data_range, options, r0set)
%
% Input parameters:
% binned_data : Binned data to be used for forming the confining potential
% data_limits : A 2xD array where the first row gives the minimum values of
%   the data and the second row gives the maximum values of the data. D is
%   the dimension of the data.
%       Ex. If the data is given by an array, dat, where each row is a
%       point and columns are dimensions, then the data_limits could be
%       computed using
%           data_limits = [min(dat,1); max(dat,1)];
% options : An instance of class seedPointOptions.
% r0set : (optional) A NxD array giving a set of possible initial
%   positions *in data space*. This is required if
%   options.Point_Selection_Method is 'r0set_random' or
%   'r0set_uniformRandom'.
%
% Output parameters:
% dV : Function handle taking in an NxD array of locations and outputting
%   the potential gradient (also an NxD array) at those locations.
% r0 : A cell array (of length options.Iterations) where each element
%   contains the initial positions to be used for the simulation.
% problem_scales : A structure with fields grid_spacing and grid_to_solver
%   giving the size of each bin in data units and the scale factor to go
%   from the grid to the solver space.
% V : The confining potential
% Info : A structure array of information that may be useful when
%   debugging. If options.Debug is false, then Info will be empty;
%   otherwise, it will have the following fields
%       ComputeInitialPointsInfo : A structure containing the debug
%           information from the computeInitialPoints() function.
%
% See also COMPUTEPROBLEMSCALES COMPUTEINITIALPOINTS

% James Kapaldo

if nargin < 4
    r0set = [];
end
Info = [];
DEBUG = options.Debug;

% Parameters
PAD_SIZE = options.Potential_Padding_Size;

% Problem scales
problem_scales = computeProblemScales(options, size(binned_data), data_limits);

% Create confining potential
switch options.Potential_Type
    case 'distance_transform'
        [V, scaleFactor] = create_scaleInvar_confining_potential(binned_data, options);
        problem_scales = scale_object_for_distance_transform(problem_scales, scaleFactor);
    case 'density'
        V = min(1.5,1./binned_data);
        V = padarray(V, PAD_SIZE*ones(1,ndims(binned_data)), 1.5);
end

% Apply potential modifier
if ~isempty(options.Potential_Modifier)
    V = options.Potential_Modifier(V);
end

% Compute confining force
dV = confiningForce(V, problem_scales, options);

% Compute initial positions
BW = (V < options.Maximum_Initial_Potential) & (V > options.Minimum_Initial_Potential);
[r0, InitPointsInfo] = computeInitialPoints(BW, options, 'problem_scales', problem_scales, 'r0set', r0set);


if DEBUG
    Info.ComputeInitialPointsInfo = InitPointsInfo;
end

end


function problem_scales = scale_object_for_distance_transform(problem_scales, scaleFactor)

dg = problem_scales.grid_spacing;
gts = problem_scales.grid_to_solver;
D = numel(dg);

% Check data and grid aspect ratios
mdg = mean(dg);
if ~all( abs(dg - mdg)/mdg < 1e-2 )
    warning('compute_ScaleFactors_and_Grid:nonEqualAspectRatios','The aspect ratio of the data and the discretized grid are not equal. This could lead to unexpected or wrong results.')
end

% We need to scale the solver space to be relative to the distance
% transform. Length on the grid (where the distance transform is defined)
% is sqrt(D), length in the data space is sqrt(dg * dg') -- *not*
% sqrt(gts*gts'). We do not use gts because gts already includes the scaling
% for the ScaleInvarient_Potential_Extent.
gts = gts .* sqrt( D / (dg*dg')); % Scale solver space so that potential extent is in units of Distance_Transform

% Now scale the grids for the object resize. Note that the distance
% transform is not related to the data limits, but only the grid; thus, the
% grid spacing, which is used for the gradient calculation, should be ones
% times the scaleFactor.
gts = gts * scaleFactor; % Scale solver space so that potential extent is in units of Max_Distance_Transform
dg = ones(1,D)*scaleFactor; % Scale grid size for gradient calculation (This effectively resizes the object so that the object's maximum distance transform value is Max_Distance_Transform.)

problem_scales.grid_spacing = dg;
problem_scales.grid_to_solver = gts;

end


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
