%% Test seed-point detection on simulated 3D data
% From the example nuclei clump used throughout the paper, create simulated
% 3D scatter point data, and then find the seed-points of the data.

%% Load data
dat = create_test_3D_object(80000);

% Scale data : the test nuclei are a ~5 times smaller in hight than diameter.
dat = dat .* [1, 1, 5];
%% Compute data density
nbins = 30;
density_threshold = 2;
smooth_data = true; % use a faster smoothing after preperation

[n, cents, sz] = binData(dat, nbins, density_threshold, smooth_data, true);

% n = frequencyGaussianFilter(n,0.5,3,'replicate');
% n(n<0) = 0;

% Get data scale factors
data_range = cellfun(@(x) range(x), cents);
data_scale = data_range./(sz-1);

%% Explore/Visualize the data

dimensions = [1,2,3];
isoLevels = [10,20,30,40];%[2,5,8,11];

ticks = cell(1,3);
ticklbls = cell(1,3);

dt = [40,40,20];
for i = 1:3
    c = cents{dimensions(i)};
    c1 = fix(c(1));
    c2 = fix(c(end));
    
    ticklbls{i} = round(c1/dt(i))*dt(i):dt(i):(c2-dt(i)/2);
    ticks{i} = interp1(c, 1:sz(dimensions(i)), ticklbls{i});
end
ticklbls{3} = ticklbls{3};
color_scale = @(x) x;

% try close(fig), catch, end;
fig = create_3d_density_plot(n, dimensions, isoLevels,'ColorScale',@(x)(x),'Ticks',ticks,'TickLabels',ticklbls);

set(gca,'XDir','reverse')
view(220,31)


%% Create mask
Omega = n > 5;

%% Locate seed-points using distance transform based confining potential
% This section should be modified. the data should be scaled better so that
% the hights are about the radii. 

options = seedPointOptions();
options.Wigner_Seitz_Radius = 2.5; % This is in units of the discretized space!
options.Wigner_Seitz_Radius_Space = 'grid';
options.Point_Selection_Method = 'uniformRandom';
options.Maximum_Initial_Potential = 1/5;
options.Potential_Padding_Size = 1;

% Scale the object so that the maximum distance transform value is 18.
options.Max_Distance_Transform = 18;

% Set the particle interaction parameter values.
options.Potential_Parameters = [-1, 2, 13];
options.Potential_Parameters_Space = 'solver'; % No

% options.Particle_Damping_Rate = 5e-4;
options.Debug = true;
options.Use_Parallel = false;

NN = 1;
seedPoints = cell(NN,1);
tmp_seedPoints = [];
for i = 1:NN
N = 20;
start = tic;

[seedPoints{i},Info] = compute_seedPoints_nd(Omega, data_range, options,'iterations',N,'minClusterSize',N/2,'verbose',0);
tmp_seedPoints = [tmp_seedPoints; Info.seedPoint_set];
fprintf('Finished, avg set time = %f\n', toc(start)/N)
end


seedPoint_set_all = cat(1,seedPoints{:});
seedPoint_set = tmp_seedPoints;%Info.seedPoint_set;
r0 = Info.iteration_info{1}.r0;

% Plot the results
options.Potential_Padding_Size = 0;
V = create_scaleInvar_confining_potential(Omega,options);
dimensions = [1,2,3];
isoLevels = [1,5,10,15];
markers = [];

markers(2).dat = Info.data_to_grid(seedPoint_set_all);
markers(2).options.Marker = '.';
markers(2).options.Color = 'r';
markers(2).options.LineStyle = 'none';
markers(2).options.MarkerSize = 15;
markers(2).project = true;

markers(1).dat = Info.data_to_grid(seedPoint_set);
markers(1).options.Marker = '.';
markers(1).options.Color = 'k';
markers(1).options.LineStyle = 'none';
markers(1).options.MarkerSize = 4;
markers(1).project = true;


try close(fig), catch, end
fig = create_3d_density_plot(1./V, dimensions, isoLevels,'ColorScale',@(x)sqrt(x), 'Markers',markers,'Ticks',ticks,'TickLabels',ticklbls);
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
%
try close(fig), catch, end
fig = create_3d_density_plot(n, dimensions, isoLevels, 'Markers', markers, 'ColorScale', @(x) (x),'Ticks',ticks,'TickLabels',ticklbls);
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