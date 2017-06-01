function [dV, r0, problem_scales, Info] = setup_problem(binned_data, data_range, options, r0set)

if nargin < 4
    r0set = [];
end
Info = [];
DEBUG = options.Debug;

% Parameters
PAD_SIZE = options.Potential_Padding_Size;

% Problem scales
problem_scales = computeProblemScales(options, size(binned_data)-1, data_range);

% Create confining potential
switch options.Potential_Type
    case 'distance_transform'
        [V, scaleFactor, overlapFactor] = create_scaleInvar_confining_potential(binned_data, options); %#ok<ASGLU>
        if ~isnan(options.Max_Distance_Transform)
            problem_scales = scale_object_for_distance_transform(problem_scales, scaleFactor);
        end
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
    Info.V = V;
    Info.ComputeInitialPoints = InitPointsInfo;
end

end


function problem_scales = scale_object_for_distance_transform(problem_scales, scaleFactor)

dg = problem_scales.grid_spacing;
gts = problem_scales.grid_to_solver;

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
dg = dg * scaleFactor; % Scale grid size for gradient calculation (This effectively resizes the object so that the object's maximum distance transform value is Max_Distance_Transform.)
%     grid_spacing = ones(1,D)*scaleFactor;

problem_scales.grid_spacing = dg;
problem_scales.grid_to_solver = gts;

end