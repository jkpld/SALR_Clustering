function problem_scales = computeProblemScales(options, grid_range, data_range)
% COMPUTEPROBLEMSCALES compute the grid spacing and the grid-to-solver
% scale factor.
%
% problem_scales = computeGridScaleFactors(options, grid_range, data_range)
%
% problem_scales : structure with two fields, grid_spacing and
% grid_to_solver.

% James Kapaldo

if nargin < 3
    data_range = [];
end

% Dimension
D = length(grid_range);

% Compute grid spacing
if isempty(data_range)
    % If not set, then assume spacing of 1
    grid_spacing = ones(1,D);
else
    grid_spacing = data_range ./ grid_range;
end

% Now scale the grid_to_solver space so that the attractive extent is
% given in terms of the ScaleInvarient_Potential_Extent
grid_to_solver = grid_spacing * options.ScaleInvarient_Potential_Extent / options.Potential_Parameters(3);

% Form output structure
problem_scales = struct('grid_spacing', grid_spacing, 'grid_to_solver', grid_to_solver);

end
