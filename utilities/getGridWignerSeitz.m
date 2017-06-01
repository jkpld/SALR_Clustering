function rs = getGridWignerSeitz(problem_scales, options)

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
