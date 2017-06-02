function [grid_spacing, grid_to_solver, grid_to_data] = computeGridScaleFactors(options, grid_range, data_range)
% COMPUTEGRIDSCALEFACTORS
%
% [grid_spacing, grid_to_solver, grid_to_data] = computeGridScaleFactors(options, grid_range, data_range)

% James Kapaldo

% The coordiante space the potential parameters are defined in
% Vint_space = options.Potential_Parameters_Space;

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


grid_to_data = grid_spacing;
grid_to_solver = grid_spacing;

% grid_spacing is the scale factors that take us from grid space to data
% space.

% If there was a scale factor defined, then scale the solver space.
if options.Scale_Factor ~= 1
    mdg = mean(grid_spacing);

    if ~all( abs(grid_spacing - mdg)/mdg < 1e-2 )
        warning('compute_ScaleFactors_and_Grid:nonEqualAspectRatios','The aspect ratio of the data and the descretized grid are not equal. This could lead to unexpected or wrong results.')
    end

    % Now we need to go from the grid to data space, but the data space
    % should be scaled so that the attractive extent is in units of the
    % Max_Distance_Transform. This is accomplished by scaling the grid
    % space to have the same aspect ratio as the data space.
    %
    % Length on the grid is sqrt(D), length in the data space is
    % sqrt(grid_spacing * grid_spacing'). Thus, the factor to take us
    % from grid to data space is
    if ~isnan(options.Max_Distance_Transform)
        grid_to_solver = grid_to_solver .* sqrt( D / (grid_spacing*grid_spacing'));
    end

    % Scale the grid spacing for the gradient calculation. The
    % scale_factor should be equal to
    %    options.Max_Distance_Transform / Object_Max_DT
    % and should have been set in just after the function that computes
    % the confining potential with the distance transfrom.

    grid_to_solver = grid_to_solver .* options.Scale_Factor; % Scale solver space so that potential extent is in units of Max_Distance_Transform
    grid_spacing = grid_spacing * options.Scale_Factor;
%     grid_spacing = ones(1,D)*options.Scale_Factor;

    % Re-set Object_Max_DT to NaN.
    options.Scale_Factor = 1;
end
% switch Vint_space
%     case 'data'
%         % Interaction potential parameters are defined in terms of the data
%
%         grid_to_solver = grid_spacing;
%
%     case 'distance_transform'
%         % Interaction potential parameters are defined in terms of the
%         % maximum value of the distance transform
%
%
%         mdg = mean(grid_spacing);
%
%         if ~all( abs(grid_spacing - mdg)/mdg < 1e-2 )
%             warning('compute_ScaleFactors_and_Grid:nonEqualAspectRatios','The aspect ratio of the data and the descretized grid are not equal. This could lead to unexpected or wrong results.')
%         end
%
%         % Now we need to go from the grid to data space, but the data space
%         % should be scaled so that the attractive extent is in units of the
%         % Max_Distance_Transform. This is accomplished by scaling the grid
%         % space to have the same aspect ratio as the data space.
%         %
%         % Length on the grid is sqrt(D), length in the data space is
%         % sqrt(grid_spacing * grid_spacing'). Thus, the factor to take us
%         % from grid to data space is
%         grid_to_solver = grid_spacing;
%         if ~isnan(options.Max_Distance_Transform)
%             grid_to_solver = grid_to_solver .* sqrt( D / (grid_spacing*grid_spacing'));
%         end
%
%         % Scale the grid spacing for the gradient calculation. The
%         % scale_factor should be equal to
%         %    options.Max_Distance_Transform / Object_Max_DT
%         % and should have been set in just after the function that computes
%         % the confining potential with the distance transfrom.
%
%         grid_to_solver = grid_to_solver .* options.Scale_Factor; % Scale solver space so that potential extent is in units of Max_Distance_Transform
%         grid_spacing = ones(1,D)*options.Scale_Factor;
%
%
%         % Re-set Object_Max_DT to NaN.
%         options.Scale_Factor = 1;
% end

% Now scale the grid_to_solver space so that the attractive extent is
% given in terms of the ScaleInvarient_Potential_Extent
grid_to_solver = grid_to_solver * options.ScaleInvarient_Potential_Extent / options.Potential_Parameters(3);








% switch scale_type
%     case 'distance_transform'
%         % This is a scale of the grid since the distance transform is
%         % defined on the grid.
%         grid_spacing = grid_spacing * scale_factor;
%         grid_to_solver_scale = scale_factor;
% %         S = DT_SCALE * grid_range ./ mean_data_range;
%     case 'data'
%         grid_to_solver_scale = grid_spacing;
% %         S = mean(d_r);
%     case 'grid'
%         grid_to_solver_scale = ones(1,D);
% %         S = grid_range ./ mean_data_range;
%     otherwise
%         grid_to_solver_scale = scale_type .* mean_data_range ./ grid_range;
% %         S = scale_type;
% end










end
