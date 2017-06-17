function rs = getGridWignerSeitz(problem_scales, options)
% GETGRIDWIGNERSEITZ Convert the Wigner-Seitz radius to grid space.
%
% rs = getGridWignerSeitz(problem_scales, options)
%
% Input parameters:
% problem_scales : The structure returned by computeProblemScales()
% options : An instance of class seedPointOptions
%
% Output parameters:
% rs : The Wigner-Seitz radius in grid space

% James Kapaldo

grid_to_data = problem_scales.grid_spacing;
grid_to_solver = problem_scales.grid_to_solver;

rs = options.Wigner_Seitz_Radius;

switch options.Wigner_Seitz_Radius_Space
    case 'data'
        % rs = rs ./ grid_to_data; % If you want a different Wigner-Seitz radius for each dimension than use this line instead of the line below.
        rs = rs / mean(grid_to_data);
    case 'grid'
    case 'solver'
        % rs = rs ./ grid_to_solver; % If you want a different Wigner-Seitz radius for each dimension than use this line instead of the line below.
        rs = rs / mean(grid_to_solver);
end

end
