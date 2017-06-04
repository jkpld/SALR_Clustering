function [dV, r0, problem_scales, V, Info] = setup_problem(binned_data, data_limits, options, r0set)
% SETUP_PROBLEM Generate the confining force, the initial particle
% locations, and the scaling factors for the simulation.
%
% [dV, r0, problem_scales, V] = setup_problem(binned_data, data_range, options)
% [dV, r0, problem_scales, V] = setup_problem(binned_data, data_range, options, r0set)
% [dV, r0, problem_scales, V, Info] = setup_problem(binned_data, data_range, options, ...)
%
% Input parameters:
% binned_data : Binned data to be used for forming the confining potential
% data_range : 1xD array giving the range of the actual data along each
%   dimension. Alternatively data_range can be empty, in which case the
%   size of each bin is assumed to be 1.
% options : An instance of class seedPointOptions.
% r0set : (optional) A NxD array giving a set of possible initial
%   positions *in data space*. This is required if
%   options.Point_Selection_Method is 'r0set_random' or
%   'r0set_uniformRandom'.
%
% Output parameters:
% dV : Function handle taking in an NxD array of locations and outputing
%   the potential gradient (also an NxD array) at those locations.
% r0 : A cell array (of length options.Iterations) where each element
%   containts the initial positions to be used for the simulation.
% problem_scales : A structure with fields grid_spacing and grid_to_solver
%   giving the size of each bin in data units and the scale factor to go
%   from the grid to the solver space.
% V : The confining potential
% Info : A structure array of information that may be useful when debuging.
%   If options.Debug is false, then Info will be empty; otherwise, it will
%   have the following fields
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
disp('problem_scales')
disp(problem_scales)
% Create confining potential
switch options.Potential_Type
    case 'distance_transform'
        [V, scaleFactor, overlapFactor] = create_scaleInvar_confining_potential(binned_data, options); %#ok<ASGLU>
        disp('scale_factor')
        disp(scaleFactor)
        if ~isnan(options.Max_Distance_Transform)
            problem_scales = scale_object_for_distance_transform(problem_scales, scaleFactor);
        end
        disp('problem_scales')
        disp(problem_scales)
    case 'density'
        V = min(1.5,1./binned_data);
        V = padarray(V, PAD_SIZE*ones(1,D), 1.5);
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
    warning('compute_ScaleFactors_and_Grid:nonEqualAspectRatios','The aspect ratio of the data and the descretized grid are not equal. This could lead to unexpected or wrong results.')
end

% We need to scale the solver space to be relative to the distance
% transform. Length on the grid (where the distance transform is defined)
% is sqrt(D), length in the data space is sqrt(dg * dg') -- *not*
% sqrt(gts*gts'). We do not use gts because gts already includes the scaling
% for the ScaleInvarient_Potential_Extent.
gts = gts .* sqrt( D / (dg*dg')); % Scale solver space so that potential extent is in units of Distance_Transform

% Now scale the grids for the object resize.
gts = gts * scaleFactor; % Scale solver space so that potential extent is in units of Max_Distance_Transform
dg = ones(1,D)*scaleFactor; % Scale grid size for gradient calculation (This effectively resizes the object so that the object's maximum distance transform value is Max_Distance_Transform.)

problem_scales.grid_spacing = dg;
problem_scales.grid_to_solver = gts;

end
