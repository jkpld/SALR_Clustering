function overlap_factor = get_overlap_factor(BW)

% Padd mask
BW = padarray(BW,[1 1]);

% Label and extract objects
L = bwlabel(BW);
CC = bwconncomp(BW);

% Get distance transform
DT = bwdist(~BW);

% Get the gradient magnitude and direction of the DT
[gdtx,gdty] = gradient(DT);
gdtm = sqrt(gdtx.^2 + gdty.^2);
gdta = atan2(gdty,gdtx);
gdta(gdta<0) = gdta(gdta<0) + pi;

% Find the local maximum of the DT and label them
m = imregionalmax(imopen(DT,strel('disk',2)));
m = bwlabel(m);

% For each local maximum, create a "ring" of radius 5 around the center
mb = imdilate(m,strel('disk',5,0));
mb = imdilate(mb,strel('disk',1)) .* (~mb);

% Get the pixels of the rings and the object numbers to which they belong
pxIdx = find(mb);
objIdx = mb(pxIdx);

% Collect all the pixels that belong to each object
pixelList = accumarray(objIdx, pxIdx, [], @(x) {x});

% Assume a line is drawn through each pixel of the rings with a slope given
% by the gradient magnitude. We will compute where this line should
% intercept DT=0 (y=0) and use the x distance to this intercept plus half
% the ring radius (plus some correction offset) to determine the radius of
% the object along the gradient direction.
radii = cellfun(@(x) DT(x)./gdtm(x) + 6+0.75*2, pixelList,'UniformOutput',false);
a = cellfun(@(x) gdta(x), pixelList,'UniformOutput',false);

% Now use the the radii to compute the area of each sub-object, modeling the
% object as an ellipse, and find the orientation of the ellipse. (A sub-
% object is what I call the object for each local maximum in the DT)
subObjAreas = cellfun(@(x,y) get_object_area(x,y), radii, a);
subObj_idx = cellfun(@(x) mean(L(x(L(x)>0))), pixelList);

% Add up all the expected sub-object areas for each object
expected_object_area = accumarray(subObj_idx(:), subObjAreas(:), [max(L(:)), 1]);
% Get the area of each object
object_area = cellfun(@numel, CC.PixelIdxList);

% Get the overlap fraction of each object
overlap_factor = expected_object_area(:) ./ object_area(:);
% overlap_factor = 1 + abs(expected_object_area(:) ./ object_area(:) - 1);
overlap_factor = double(overlap_factor);

end


function area = get_object_area(radii,a)

% Bin the angles
a_edge = linspace(0,pi+0.1,10);
bin = discretize(a,a_edge);

% Average the radii over each angular bin
radii = accumarray(bin(all(~isnan(bin),2)), radii(:), [numel(a_edge)-1,1], @mean);
radii(radii==0) = [];

% sort the radii and assume an ellipse with minor radius given by the
% smallest radii and the major radius given by the 75th percentile radii.
sr = sort(radii);
area = pi * sr(1) * interp1(linspace(0,1,numel(sr)), sr, 0.75);

end
