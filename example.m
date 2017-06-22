%% Example script for declumping objects
% Make sure to first run setup.m


% Get a test image =======================================================
I = imread('testImage_image_LD4P24.tif');
BW = imread('testImage_mask_LD4P24.tif');

% Intialize options ======================================================
options = declumpOptions();

options.Max_Radius = 35;
options.Min_Angle = 0.5;

options.Wigner_Seitz_Radius = 20;
options.Potential_Depth = -1;
options.Potential_Minimum_Location = 2;
options.Potential_Extent = 15;

options.Point_Selection_Method = 'curvatureUniformRandom';

options.Use_GPU = false;
options.Use_Parallel = false;

options.Debug = true; % Set this to false if you only need the declumped mask. Leave as true to get the extra information.

    % Note: if you want to see details of a particular object then
    % uncomment the line below and put in the object number you are
    % interested in.
    % options.Object_Of_Interest = 1;

% Declump the objects in the image =======================================
[declumpedBW, cuts, Info] = declumpNuclei(I,BW,options);

% Note: I have not yet writen up a help section for the declumpNuclei
% function yet, but it is commented. declumpNuclei calls the primary
% function declumpObject__, and the functions that this function calls do
% have the help sections.


% Show the segmentation as rgb image =====================================
CC = bwconncomp(declumpedBW,8);
L = labelmatrix(CC);
RGB = label2rgb(L,'jet',[0 0 0],'shuffle');
fig = figure('visible','off');
imshow(RGB)
hold on
ax = gca;
title('declumped objects')
drawnow;
goDark(gcf);

% Add the nuclei centers to the image as white dots
th = linspace(0,360,50);
r = 3;
x = r*cosd(th);
y = r*sind(th);
patchCircle = @(center) patch('XData',x+center(2),'YData',y+center(1),'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor',[1 1 1],'Parent',ax,'HitTest','off');
if ~isempty(Info)
    nonEmptyInfo = cellfun(@(x) isfield(x,'centers'),Info);
    centers = cellfun(@(x) x.centers, Info(nonEmptyInfo),'UniformOutput',false);
    centers = cat(1,centers{:});
    for i = 1:size(centers,1)
        patchCircle(centers(i,:))
    end
%     line(centers(:,2),centers(:,1),'Marker','o','LineStyle','none','Color',[0 0 0],'MarkerSize',3,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0 0 0])
end
fig.Visible = 'on';