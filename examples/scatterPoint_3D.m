%% Test seed-point detection on simulated 3D data
% From the example nuclei clump used throughout the paper, create simulated
% 3D scatter point data, and then find the seed-points of the data.

%% Load data
dat = create_test_3D_object(80000);

% Scale data: In order to use the distance transform each object should be
% about the size size along any dimension. The objects in this data set are
% about 5 times smaller in height than their latteral size; so, scale the
% heights by 5.
dat = dat .* [1, 1, 5];


%% Bin the data
nbins = 30;
density_threshold = 2;
smooth_data = true;
[n, cents] = binData(dat, nbins, density_threshold, smooth_data, true);

%% Explore/Visualize the data
dimensions = [1,2,3];
isoLevels = [10,20,30,40];
tick_spacing = [30,30,20];

fig = create_3d_density_plot(n, dimensions, isoLevels, ...
    'bin_centers', cents, 'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)

%% Create mask
Omega = n > 5;

%% Locate seed-points using distance transform based confining potential
options = seedPointOptions();

% Set the particle initialization parameters
options.Point_Selection_Method = 'uniformRandom';
options.Wigner_Seitz_Radius = 2.5;
options.Wigner_Seitz_Radius_Space = 'grid';
options.Maximum_Initial_Potential = 1/5;

% Set the particle interaction parameter values.
options.Potential_Parameters = [-1, 2, 13];

% Scale the object so that the maximum distance transform value is 18.
options.Max_Distance_Transform = 18;
options.Potential_Padding_Size = 2;

% Set up replicates and minimum cluster size.
options.Iterations = 10;
options.Minimum_Cluster_Size = 10/3;

options.Verbose = true;
options.Debug = true;
options.Use_Parallel = false;

% Repeat process NN times to check reproducability.
NN = 1;
seedPoints = cell(NN,1);
seedPoints_n = [];

for i = 1:NN
    start = tic;

    % Compute the seed points
    [seedPoints{i},Info] = computeObjectSeedPoints(Omega, options,'data_limits',data_limits);

    % Save the seed point results of each replicate
    seedPoints_n = [seedPoints_n; Info.seedPoints_n]; %#ok<AGROW>

    fprintf('Finished, avg set time = %f\n', toc(start)/options.Iterations)
end

seedPoints = cat(1,seedPoints{:});

% Plot the results
V = removePadding(Info.V, options.Potential_Padding_Size);

dimensions = [1,2,3];
isoLevels = [1,5,10,15];
markers = [];

% Helper function to convert are seed points into grid space
data_to_grid = @(t) Info.problem_scales.data_to_grid(t) - ...
    options.Potential_Padding_Size;

% Plot the final seed points as large red dots
markers(2).dat = data_to_grid(seedPoints);
markers(2).options.Marker = '.';
markers(2).options.Color = 'r';
markers(2).options.LineStyle = 'none';
markers(2).options.MarkerSize = 15;
markers(2).project = true;

% Plot the seed points of each replicate as small black dots
markers(1).dat = data_to_grid(seedPoints_n);
markers(1).options.Marker = '.';
markers(1).options.Color = 'k';
markers(1).options.LineStyle = 'none';
markers(1).options.MarkerSize = 4;
markers(1).project = true;

try close(fig), catch, end
fig = create_3d_density_plot(1./V, dimensions, isoLevels, ...
    'ColorScale', @(x)sqrt(x), ...
    'Markers', markers, ...
    'bin_centers', cents, ...
    'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view(220,31)

%% Locate seed-points with k-means
rng('default') % ensure we get the same kmeans results each time
N = 20;
c = cell(N,1);

kmeans_options = [];
kmeans_options.MaxIter = 200;
kmeans_options.UseParallel = false;
tic
for i = 1:N
    [~,c{i}] = kmeans(dat,7,'Start','plus','Options',kmeans_options,'Replicates',2);
end
fprintf('time per kmeans run : %f\n', toc/N)
c = cat(1,c{:});

for i = 1:3
    c(:,i) = interp1(cents{i},1:size(n,i),c(:,i));
end

% Plot the results
dimensions = [1,2,3];
isoLevels = 10:10:40;

markers = [];
markers(1).dat = c;
markers(1).options.Marker = '.';
markers(1).options.Color = 'r';
markers(1).options.LineStyle = 'none';
markers(1).options.MarkerSize = 20;
markers(1).project = false;

try close(fig), catch, end
fig = create_3d_density_plot(n, dimensions, isoLevels, ...
    'Markers', markers, ...
    'bin_centers', cents, ...
    'tick_spacing', tick_spacing);

set(gca,'XDir','reverse')
view([180,90])
%% Export figure
% Export lines as pdf.
% Export transparent isosurfaces as png.
% Overlay the png on the pdf in other software.
%
% Note, this is distructive! It will delete the transparent objects before
% the end.
%
% The code uses export_fig for exporting the images.
% https://github.com/altmany/export_fig
% https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

fig.Units = 'inch';
fig.PaperSize = fig.Position(3:4);
fig.PaperPosition = [0 0 fig.Position(3:4)];
fig.PaperUnits = 'inch';

pth = 'K:\Google_Drive\NotreDame\RadiationLab\Projects\nucleiDeclumping\Manuscript\Figures\scatterPoint_figure\results_20_iterations';

h = [findall(fig,'Type','line'); findall(fig,'UserData','contourLine')];
htxt = findall(fig,'Type','text');
hm = findall(fig,'UserData','isoSurfaceMarker'); % This is a marker to help properly overlay the images.

set([h;htxt],'Visible','off')
set(hm,'Visible','on')

% Export the transparent isosurfaces
% export_fig(fig,'-png','-r400','-transparent',pth)

set([h;htxt],'Visible','on')

hiso = findall(fig,'Tag','isoSurface');
delete(hiso)
delete(hm)

% Export the lines
% export_fig(fig,'-painters','-pdf',pth)

% Export xy plane
view([180,90])
export_fig(fig,'-painters','-pdf',[pth, '_xy'])
