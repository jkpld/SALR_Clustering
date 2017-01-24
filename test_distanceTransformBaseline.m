%% Base line : distance transform -> regionalmax -> mean-shift
% Compare results of new method with base line distance transform. Since
% the new method uses the distance transform, I want to see if the new
% method is good, or if the distance transform is good.

pth = 'K:\Google_Drive\MATLAB\seed_point_detection';

im_pth      = @(n) [pth '\exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth '\exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth '\exampleImages\markedCenters_LD' n 'P24'];

names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10;

TP = zeros(n,numel(dr));
TPplusFP = zeros(n,1);
TPplusFN = zeros(n,1);
dN = zeros(n,1,7);

for ind = 1:n
    fprintf('%s >> Image = LD%sP24 (%d/%d)...\n', datestr(now,31),names{ind},ind,numel(names))
    BW = imread(bw_pth(names{ind})) > 0;
    L = bwlabel(BW);

    results = load(results_pth(names{ind}));
    result_fields = fieldnames(results);
    tmp = results.(result_fields{1});

    NS = createns(tmp(:,[3,2]));
    correctNumCents = accumarray(tmp(:,1),1,[],[],0);
    objNumbers = unique(tmp(:,1));
    trueNumberOfNuclei = size(tmp,1);
    TPplusFN(ind,:) = trueNumberOfNuclei;

    D = bwdist(~BW);
    D = imopen(D,strel('disk',2));
    M = imregionalmax(D);
    props = regionprops(M,'Centroid');
    seedPoints = cat(1,props.Centroid);
    objNum = round(interp2mex(L,seedPoints(:,1),seedPoints(:,2)));

    for drNo = 1:numel(dr)
        idx = rangesearch(NS,seedPoints,dr(drNo));
        TP(ind,drNo) = numel(unique(cat(2,idx{:})));
    end

    TPplusFP(ind) = size(seedPoints,1);

    %Number of centers difference
    numCents = accumarray(objNum,1,[max(objNumbers),1],[],0);
    d = numCents(objNumbers) - correctNumCents(objNumbers);

    d(abs(d)>3) = [];
    d = d+ 4;
    dN(ind,:) = accumarray(d,1,[7,1],[],0);
end

Pnn = TP./TPplusFP;
Rnn = TP./TPplusFN;
F1nn = 2*Pnn.*Rnn./(Pnn+Rnn);

% Results for all images
P = sum(TP,1)./sum(TPplusFP,1);
R = sum(TP,1)./sum(TPplusFN,1);
F1 = 2*P.*R./(P+R);

results = [];
results.Pnn = Pnn;
results.Rnn = Rnn;
results.F1nn = F1nn;
results.P = P;
results.R = R;
results.F1 = F1;
results.dN = dN;
results.TP = TP;
results.TPplusFP = TPplusFP;
results.TPplusFN = TPplusFN;

dims = [];
dims.PnnRnnF1nn = {'image','dr'};
dims.PRF1 = {'1','dr'};
dims.dN = {'image','dNm3top3'};
dims.TP = {'image','dr'};
dims.TPplusFPorFN = {'image'};
dims.dr = dr;
dims.dN = -3:3;
dims.image = {'LD2P24','LD3P24','LD4P24','LD5P24','LD67P24'};
results.dims = dims;
results.note = 'bwdist -> imopen(disk,2) -> imregionalmax -> centroid';


fprintf('Distance transform baseline, open(disk,2)\n')
disp(squeeze(F1))
disp(squeeze(sum(dN,1))'/(5*484))