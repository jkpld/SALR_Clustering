%% 3D scatter-point data
% This example uses a simple 3D data set, made to resemble 3D nuclei, for
% the purpose of validating the SALR particle clustering result with the
% k-means clustering result. The confining potential will be based on the
% distance transform, and all parameters for the SALR clustering, other
% than the particle's initial locations, will be exactly the same as those
% used for the 2D nuclei example.

%% Create data set
dat = create_test_3D_object(80000);

%% Bin the data
% In this step the data is discretized to a grid (binned). The data is
% initially binned so that the grid completely covers the data range. Using
% this information the grid range is limited so that the bins have at least
% `density_threshold` points. The data is then re-binned using this limited
% data range and smoothed with a slowish, but memory efficient Gaussian
% filter. The grid is set to have approximately `nbins` along each
% dimension. The actual number of bins will be chosen so that the grid has
% the same aspect ratio as the data.

nbins = 30;
density_threshold = 2;
smooth_data = true;
[n, cents] = binData(dat, nbins, density_threshold, smooth_data, true);

%% Explore/Visualize the data
dimensions = [1,2,3];
isoLevels = [10,20,30,40];
tick_spacing = [40,40,8];

fig = isosurfaceProjectionPlot(n, dimensions, isoLevels, ...
    'bin_centers', cents, 'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)
% <INSERT FIGURE> -- this tag is used in generating the published output
% <CAPTION>
% Note the scale difference between the object heights and lateral sizes.
% </CAPTION>

%%
% The figure shows that the height of the objects is about 5 times smaller
% than the lateral sizes. If the data scale was left as it is, the
% distance transform would not work well as there would be very little
% force in the lateral directions. Thus, the data should first be scaled.
%
%#### Scale data and re-bin

dat = dat .* [1, 1, 5];
[n, cents, sz, data_limits] = binData(dat, nbins, density_threshold, smooth_data, true);
tick_spacing = [40,40,20];

fig = isosurfaceProjectionPlot(n, dimensions, isoLevels, ...
    'bin_centers', cents, 'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)
% <INSERT FIGURE> -- this tag is used in generating the published output
% <CAPTION>
% Note the scale of the objects in all three dimensions is now about the
% same.
% </CAPTION>

%% Create mask
% To base the confining potential on the distance transform, a binary mask
% is needed. Here the data density is thresholded at 5 to create the mask.
Omega = n > 5;

%% Setup SALR clustering parameters
% The parameters needed to run SALR clustering can all be accessed and set
% using the class `seedPointOptions`. This class handles parameter
% validation as well as computing/visualizing the SALR particle interaction
% parameters/potential.
options = seedPointOptions();

%%
%#### Set the particle initialization parameters.
%
% * Use a uniform random point distribution. This will overlay a
% hyper-cubic lattice across the grid with each lattice cell having a
% volume equal to the volume of a hyper-sphere with radius
% `Wigner_Seitz_Radius`. From each lattice cell, a point is then randomly
% selected from the grid where the confining potential is less than
% `Minimum_Initial_Potential`.
options.Point_Selection_Method = 'uniformRandom';
options.Wigner_Seitz_Radius = 2.5;
options.Wigner_Seitz_Radius_Space = 'grid';
options.Maximum_Initial_Potential = 1/5;

%%
%#### Set confining potential parameters.
%
% * Use a confining potential based on the distance transform
% * Scale the object so that its maximum distance transform value is 18.
% * Add a small padding around the object to ensure the confining potential
% defined around the object.
options.Potential_Type = 'distance_transform';
options.Max_Distance_Transform = 18;
options.Potential_Padding_Size = 2;

%%
%#### Set the particle interaction parameter values.
%
% * _Note_ the `Potential_Parameters` are given in data units after scaling
% the data space so that the object has a maximum distance transform value
% of `Max_Distance_Transform`.
options.Potential_Parameters = [-1, 2, 13];

%%
%#### Set up replicates and minimum cluster size.
%
% * Here we use 10 replicates and we keep any seed-point that at least 4 of
%   the replicates produce.
options.Iterations = 10;
options.Minimum_Cluster_Size = 10/3;

%%
%#### Set parameters controlling execution.
%
% * Verbose will output information on the current iteration and the
% expected time of completion.
% * Debug will return extra information.
% * Use_Parallel will determine if the iterations are run in parallel.
options.Verbose = true;
options.Debug = true;
options.Use_Parallel = false;

%% Compute seed points
% The seed-points are simply computed by passing the binned data, the
% `seedPointOptions`, and the data limits.
[seedPoints, Info] = computeObjectSeedPoints(Omega, options, 'data_limits', data_limits);

%% Plot the results
% Plot the the final seed-points as large red dots and the seed-points from
% each repetition as small black dots. This can be done by creating the
% markers structure below and passing it to the `isosurfaceProjectionPlot`
% function. It is also nice to project the seed-points onto the three axis
% planes; this can be done by setting a `project` field in the markers
% structure to true.
% * If `options.Debug = true`, then the confining potential will be
% returned in the `Info` structure. This will be used when plotting the
% results instead of the data density.
% * _Note_: The seed-points returned by `computeObjectSeedPoints` are in
% _data units_. In order to plot them, we need to convert them back into
% _grid units_.

% Remove the padding from the confining potential
V = removePadding(Info.V, options.Potential_Padding_Size);

dimensions = [1,2,3];
isoLevels = [1,5,10,15];
markers = [];

% Helper function to convert are seed points into grid space
data_to_grid = @(t) Info.problem_scales.data_to_grid(t) - ...
    options.Potential_Padding_Size;

lineSpec = @(col,ls,m,ms) struct('Color', col, ...
'LineStyle', ls, 'Marker', m, 'MarkerSize', ms);

markers(3).dat = data_to_grid(seedPoints);
markers(3).options = lineSpec('r','none','.',15);
markers(3).project = true;

markers(1).dat = data_to_grid(Info.seedPoints_n);
markers(1).options = lineSpec('k','none','.',4);
markers(1).project = true;

fig = isosurfaceProjectionPlot(1./V, dimensions, isoLevels, ...
    'ColorScale', @(x)sqrt(x), ...
    'Markers', markers, ...
    'bin_centers', cents, ...
    'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)
% <INSERT FIGURE> -- this tag is used in generating the published output
% <CAPTION>
% Small black markers represent the seed-points calculated by each
% repetition. Large red markers represent the final seed-points.
% </CAPTION>

%% Locate seed-points with [k-means][1]
% Since this data set is quite simple, we expect [k-means][1] clustering to
% work well. Thus, we can use it to validate the SALR clustering results.
% * If the number of clusters used with [k-means][1] is changed, it can be
% found that using `K=7` clusters is the most stable (see the
% publication[^1]). That will be the number of clusters used here, and we
% will use two replicates.
% * The [k-means][1] and SALR clustering results will be compared by
% plotting the [k-means][1] results as large blue dots with the SALR
% clustering results.

K = 7;
kmeans_options = [];
kmeans_options.MaxIter = 200;
kmeans_options.UseParallel = false;

[~,c] = kmeans(dat,K,'Start','plus',...
                     'Options',kmeans_options,...
                     'Replicates',2);

% Plot the results and compare with SALR particle clustering
markers(2).dat = data_to_grid(c);
markers(2).options = lineSpec('b','none','.',15);
markers(2).project = true;

fig = isosurfaceProjectionPlot(1./V, dimensions, isoLevels, ...
    'ColorScale', @(x)sqrt(x), ...
    'Markers', markers, ...
    'bin_centers', cents, ...
    'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)
% <INSERT FIGURE, THUMB> -- this tag is used in generating the published output
% <CAPTION>
% Large red markers represent the final seed-points of SALR particle
% clustering. Large blue markers represent the results of [k-means][1]
% clustering. <b>SALR clustering is able to reproduce k-means
% clustering.</b>
% </CAPTION>

%%
% *[SALR]: short-range attractive long-range repulsive
% [^1]: J. Kapaldo et al. (submitted)
% [1]: https://en.wikipedia.org/wiki/K-means_clustering
