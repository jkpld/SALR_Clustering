function dat = create_test_3D_object(N)
% CREATE_TEST_3D_OBJECT Create a fake set of 3D data for testing based off
% of the example nuclei image, testNuclei_mask.
%
% dat = create_test_3D_object(N)
%
% Input parameters:
% N : number of points to generate
%
% Output parameters:
% dat : Nx3 array of points

% James Kapaldo

if nargin < 1
    N = 20000;
end

% Read in image ----------------------------------------------------------

pth = 'exampleData\';
I = imread([pth, 'testNuclei_image.tif']);
BW = imread([pth, 'testNuclei_mask.tif']);
declumpedBW = imread([pth, 'testNuclei_segmentedMask.tif']);

BWd = imdilate(BW,strel('disk',5));
I = double(I);
I = I .* double(BWd);


% Add height -------------------------------------------------------------

Iz = zeros([size(I),15]);
Iz(:,:,8) = I;
erosionSize = [1,1,2,4,5,5,5];
tmpI = I;
for i = 1:7
    tmpI = 0.8*imerode(tmpI,strel('disk',erosionSize(i)));
    Iz(:,:,8-i) = tmpI;
    Iz(:,:,8+i) = tmpI;
end

% Create probability distribution from 3d intensity ---------------------
P = double(Iz)/trapz(trapz(trapz(Iz)));
P = P(:);
CP = cumsum(P);

[X,Y,Z] = meshgrid(1:size(I,2),1:size(I,1),1:size(Iz,3));

% offset the nuclei so they are not on the same z-level
nucNum = bwlabel(declumpedBW);
L = zeros(size(nucNum));
offset = [2,0,-2,-2,2,0,2,-2];

for i = 1:max(nucNum(:))
    L = L + offset(i)*(nucNum == i);
end

L = imfilter(L,fspecial('gaussian',21,3));

Z = Z + L;

X = X(:);
Y = Y(:);
Z = Z(:);

% Take N numbers from the distribution -----------------------------------
rng('default') % make sure we use the same random numbers for testing each time.
nP = rand(N,1);
try
    nIdx = nakeinterp1(CP, (1:numel(CP))', nP);
catch
    nIdx = interp1(CP, (1:numel(CP))', nP);
end
nIdx = round(nIdx);

dat = [X(nIdx),Y(nIdx),Z(nIdx)];
dat = dat + 0.2*randn(N,3).*std(dat);

end