%% ======================================================================
pth = 'exampleImages\';
I = imread([pth, 'testNuclei_image.tif']);
BW = imread([pth, 'testNuclei_mask.tif']);
BW = imfill(BW,'holes');

I(1:2,:) = [];
BW(1:2,:) = [];

I(end-1:end,:) = [];
BW(end-1:end,:) = [];
I(:,end-2:end) = [];
BW(:,end-2:end) = [];
BW = logical(BW);

options = seedPointOptions();
options.Wigner_Seitz_Radius = 5;
options.Debug = true;

[seedPoints,Info] = computeNucleiCenters([],BW,options);
figure(10)
clf(10)
imshow(BW)
hold on
plot(seedPoints(:,2),seedPoints(:,1),'.r')

Info{1}