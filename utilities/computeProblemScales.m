function problem_scales = computeProblemScales(options, grid_size, data_limits)
% COMPUTEPROBLEMSCALES Compute the scale factors and data transformations
% for the problem.
%
% problem_scales = computeProblemScales(options, grid_size, data_limits)
%
% Input parameters:
% options : An instance of class seedPointOptions.
% grid_size : The size of the grid.
% data_limits : A 2xD array where the first row gives the minimum values of
%   the data and the second row gives the maximum values of the data. D is
%   the dimension of the data.
%       Ex. If the data is given by an array, dat, where each row is a
%       point and columns are dimensions, then the data_limits could be
%       computed using
%           data_limits = [min(dat,1); max(dat,1)];
%
% Output parameters:
% problem_scales : structure with four fields
%   grid_spacing : 1xD array giving the size of each grid bin in terms of
%       the data. D is the dimension of the data.
%   data_to_grid : A function handle that converts from data space to grid
%       space, while also taking the Potential_Padding into account.
%   grid_to_data : A function handle that is the inverse transform of
%       data_to_grid.
%   grid_to_solver : 1xD array giving the conversion factors to go from
%       grid space to solver space.
%
% Notes: 
% * If data_limits is not given or is empty, then data_limits will
%   be set to data_limits = [zeros(1,D), grid_size]. This will give a
%   grid_spacing of 1.
%
% * If options is empty, then the Potential_Padding_Size will be
%   considered to be 0, and the attractive extents will be considered to be
%   1.
%
% See also SEEDPOINTOPTIONS

% James Kapaldo

% Options.
if isempty(options)
    PAD_SIZE = 0;
    solver_ra = 1;
    data_ra = 1;
else
    PAD_SIZE = options.Potential_Padding_Size;
    solver_ra = options.ScaleInvarient_Potential_Extent;
    data_ra = options.Potential_Parameters(3);
end

% Dimension
D = length(grid_size);

if (nargin < 3) || isempty(data_limits)
    % If the data limits are not given, assume the grid spacing is 1 and
    % there is no translation between grid and data spaces. Thus, the data
    % limits are
    data_limits = [zeros(1,D); grid_size];
end

% Compute grid spacing
%   The grid spacing is used when computing the gradient of the confining
%   potential.
data_range = diff(data_limits,1,1);
grid_spacing = data_range ./ grid_size;

% Compute grid <--> data conversions.

data_to_grid = @(r) (r - data_limits(1,:))./grid_spacing + PAD_SIZE;
grid_to_data = @(r) (r - PAD_SIZE).*grid_spacing + data_limits(1,:);

% Compute grid --> solver conversion.
%   Note that solver space is simply the data space scaled by a factor. The
%   factor should be such that the attractive extent in solver space is
%   given by ScaleInvarient_Potential_Extent
grid_to_solver = grid_spacing * solver_ra / data_ra;

% Form output structure
problem_scales = struct('grid_spacing', grid_spacing, ...
                        'data_to_grid', data_to_grid, ...
                        'grid_to_data', grid_to_data, ...
                        'grid_to_solver', grid_to_solver);

end


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
