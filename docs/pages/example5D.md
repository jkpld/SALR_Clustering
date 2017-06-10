---
layout: code-example
title:  "Example with 5D data"
subheadline:  "SALR clustering"
teaser: "The data used in this example are features representing the amount of damage in the nuclei used throughout this work. The first feature is the log of the amount of DNA double strand break normalized by the about of DNA for each nuclei. The other four features are the first four principle components (PCA) of texture and granularity features from the double strand break image channel (the images giving the double strand breaks are not provided with the example data of this work.)"
tags:
    - design
    - background color
    - header
permalink: /example5d/
---

## Test seed-point detection on real 5D data
The data used in this example are features representing the amount of
damage in the nuclei used throughout this work. The first feature is the
log of the amount of DNA double strand break normalized by the about of
DNA for each nuclei. The other four features are the first four
principle components (PCA) of texture and granularity features from the
double strand break image channel (the images giving the double strand
breaks are not provided with the example data of this work.)
## Load data

{% highlight matlab %}
dat = load('exampleImages\damage_data_5D.mat');
dat = dat.dat;

% Scale by stddev
dat_std = std(dat);
dat = dat./dat_std;
{% endhighlight %}

## Compute data density

{% highlight matlab %}
nbins = 30;
density_threshold = 20;
smooth_data = true; % Slowish but memory efficient smoothing with Gaussian: sigma=0.5, kernal_size=3

[n, cents, sz, data_limits] = binData(dat, nbins, density_threshold, smooth_data, 1);
{% endhighlight %}

## Explore/Visualize the data

{% highlight matlab %}
dimensions = [1,3,4];
isoLevels = [20,150,500,2000];
tick_spacing = 2*ones(1,5);

cmap = brewermap(9,'GnBu');
cmap(1:2,:) = [];

create_3d_density_plot(n, dimensions, isoLevels,...
'ColorScale',@(x) x.^(1/2), ...
'Bin_Centers', cents, ...
'Tick_Spacing', tick_spacing, ...
'ColorMap',cmap);

view([50,35])
% daspect(1./data_range(dimensions))
{% endhighlight %}

## Setup SALR clustering parameters

{% highlight matlab %}
options = seedPointOptions();

% Set the particle initialization parameters
options.Point_Selection_Method = 'uniformRandom';
options.Wigner_Seitz_Radius = 5;
options.Wigner_Seitz_Radius_Space = 'grid';
options.Maximum_Initial_Potential = 1/4;
options.Minimum_Initial_Potential = 1/6;

% Set confining potential parameters
options.Potential_Type = 'density';
options.Max_Potential_Force = 0.4;
options.Potential_Padding_Size = 0;
options.Maximum_Memory = 2; % Allow 2 GB for potential gradients per worker.

% Set the particle interaction parameter values.
options.Solver_Space_Attractive_Extent = 15;
options.Potential_Parameters = [-1, 0.15*2, 2];
options.Distance_Metric = {'min',4}; % Minkowski distance metric with exponent 4

% Set up replicates and minimum cluster size.
options.Iterations = 5;
options.Minimum_Cluster_Size = 3;%options.Iterations/3;

options.Verbose = true;
options.Debug = true;
options.Use_Parallel = false;

{% endhighlight %}

## Compute seed points

{% highlight matlab %}
rng('shuffle') % Ensure in random state
[seedPoints,Info] = computeObjectSeedPoints(n, options, 'data_limits', data_limits);

{% endhighlight %}

## Plot the results

{% highlight matlab %}
dimensions = [1,3,4];
isoLevels = [20,150,500,2000];
markers = [];

% Helper function to convert are seed points into grid space
data_to_grid = @(t) Info.problem_scales.data_to_grid(t) - ...
options.Potential_Padding_Size;

markers(3).dat = data_to_grid(Info.seedPoints_n);
markers(3).options.Marker = '.';
markers(3).options.Color = 'k';
markers(3).options.LineStyle = 'none';
markers(3).options.MarkerSize = 4;
markers(3).project = true;

markers(2).dat = data_to_grid(seedPoints);
markers(2).options.Marker = '.';
markers(2).options.Color = 'r';
markers(2).options.LineStyle = 'none';
markers(2).options.MarkerSize = 15;
markers(2).project = true;

% markers(1).dat = Info.iteration_info{1}.r0;
% markers(1).options.Marker = '.';
% markers(1).options.Color = 'k';
% markers(1).options.LineStyle = 'none';
% markers(1).options.MarkerSize = 4;
% markers(1).project = false;
%
try close(fig), catch, end
create_3d_density_plot(n, dimensions, isoLevels,...
'ColorScale',@(x) x.^(1/2), ...
'Markers', markers, ...
'Bin_Centers', cents, ...
'Tick_Spacing', tick_spacing, ...
'ColorMap',cmap);

view([50,35])
{% endhighlight %}

## Locate seed-points with k-means

{% highlight matlab %}
% rng('default')
% N = 1;
% c = cell(N,1);
%
% kmeans_options = [];
% kmeans_options.MaxIter = 200;
% kmeans_options.UseParallel = false;
%
% kmeans_times = zeros(1,N);
% for i = 1:N
%     tic
%     [~,c{i}] = kmeans(dat,5,'Start','plus','Options',kmeans_options,'Display','final','Replicates',2);
%     kmeans_times(i) = toc;
%     kmeans_times(i)
% end
% [mean(kmeans_times), std(kmeans_times)]
% c = cat(1,c{:});
%
% for i = 1:5
%     c(:,i) = interp1(cents{i},1:size(n,i),c(:,i));
% end
%
% % Plot the results
%
% dimensions = [1,3,4];
% isoLevels = [20,150,500,2000];
%
% markers = [];
%
% markers(1).dat = c;
% markers(1).options.Marker = '.';
% markers(1).options.Color = 'r';
% markers(1).options.LineStyle = 'none';
% markers(1).options.MarkerSize = 12;
% markers(1).project = true;
%
% try close(fig), catch, end
% fig = create_3d_density_plot(n, dimensions, isoLevels, 'Markers', markers, 'ColorScale', @(x) sqrt(x),'Ticks',ticks,'TickLabels',ticklbls,'Colormap',cmap);
% view([50,35])
{% endhighlight %}

## Export figure
