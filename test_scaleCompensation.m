%% Test restructured code
pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];


% Initialize options ----------------------------------------------------

WIGNER_SEITZ_RADIUS = 5;
CURVATURE_MAX_RADIUS = 35;
CURVATURE_SMOOTHING_SIZE = 2;
INITIAL_SPEED = 0.01;
MINIMUM_HOLE_SIZE = 15;

options = seedPointOptions();
options.Wigner_Seitz_Radius         = WIGNER_SEITZ_RADIUS;
options.Potential_Depth             = -1;
options.Potential_Minimum_Location  = 2;
options.Potential_Extent            = 15;
options.Minimum_Hole_Size           = MINIMUM_HOLE_SIZE;
options.Use_GPU                     = true;
options.Use_Parallel                = true;
options.Point_Selection_Method      = 'r0set_uniformRandom';
options.Debug                       = true;

options.Potential_Scale             = 19;
% options.Maximum_Initial_Potential   = 1/2;
% options.Object_Of_Interest          = 1;

% options.Particle_Damping_Rate = 5e-2;

% Set up test parameters ------------------------------------------------
names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10;
N = 1;

% [X,Y] = meshgrid(0.7:0.15:2, 6:2:30);
% Y(:,2:2:end) = Y(:,2:2:end) + 1;
% ra = Y(:);
% s = X(:);
% 
% r0 = round(ra*0.12);
% % ra = 5:5:30;
% % ps = 5:5:35;

[X,Y] = meshgrid(0.7, 6:4:34);
Y(:,1:2:end) = Y(:,1:2:end) + 2;
ra = Y(:);
s = X(:);


r0 = round(ra*0.12);

% [ra,ps] = meshgrid(ra,ps);
% params = [ra(:),s(:),r0(:)];
params = [14,1,2];
options.Potential_Scale = 16;
% Initialize variables for storing resuts -------------------------------
TP = zeros(n,N,size(params,1),numel(dr));
TPplusFP = zeros(n,N,size(params,1));
TPplusFN = zeros(n,N,size(params,1));
dN = zeros(n,N,size(params,1),7);

Info = struct('centers',[],'solverTime',[],'totalComputationTime',NaN,'N',NaN,'message',[]);
Info(n,N,size(params,1)).totalComputationTime = NaN;
% Info = cell(484,1);

% Compute results -------------------------------------------------------

for ind = 1:n

    fprintf('%s >> Image = LD%sP24 (%d/%d)...\n', datestr(now,31),names{ind},ind,numel(names))
    
    I = imread(im_pth(names{ind}));
    BW = imread(bw_pth(names{ind})) > 0;
    L = bwlabel(imfill(I~=0,'holes'));
    
    results = load(results_pth(names{ind}));
    result_fields = fieldnames(results);
    tmp = results.(result_fields{1});
    truthDat = tmp;
    
    
    correctNumCents = accumarray(tmp(:,1),1,[],[],0);
    objNumbers = unique(tmp(:,1));
    trueNumberOfNuclei = size(tmp,1);
    TPplusFN(ind,:,:) = trueNumberOfNuclei;

    for pi = 1:size(params,1)
        if ~mod(pi,2)
            fprintf('  %s >> param %d/%d...\n', datestr(now,31),pi,size(params,1))
        end
        
        if (strcmp(names{ind},'67') || strcmp(names{ind},'5')) && params(pi,2)>4
            options.Use_Parallel = false;
        else
            options.Use_Parallel = true;
        end
        
        options.Potential_Minimum_Location = params(pi,3);
        options.Potential_Extent = params(pi,1);
        
        options.Wigner_Seitz_Radius = round( params(pi,2) * WIGNER_SEITZ_RADIUS);
        options.Curvature_Max_Radius = ceil( params(pi,2) * CURVATURE_MAX_RADIUS);
        options.Curvature_Smoothing_Size = max(1, ceil( params(pi,2) * CURVATURE_SMOOTHING_SIZE));
        options.Initial_Speed = params(pi,2) * INITIAL_SPEED;
        
        options.Minimum_Hole_Size = round( params(pi,2) * MINIMUM_HOLE_SIZE);

        BWs = imresize(BW,params(pi,2),'bilinear','Antialiasing',false)>0.5;
        NS = createns(truthDat(:,[3,2])*params(pi,2));
        drs = dr*params(pi,2);
        
        for ni = 1:N
            if ~mod(ni,10)
                fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
            end
    %         profile on
            start = tic;
            [seedPoints, runInfo] = computeNucleiCenters_distTransform(I,BWs,options);
            totalTime = toc(start);
    %         profile off
            Info(ind,ni,pi).centers = seedPoints;
            Info(ind,ni,pi).solverTime = cellfun(@(x) x.solverTime, runInfo);
            Info(ind,ni,pi).totalComputationTime = totalTime;
            Info(ind,ni,pi).N = cellfun(@(x) size(x.r0,1), runInfo);
            Info(ind,ni,pi).message = cellfun(@(x) x.message, runInfo);

%             objNum = seedPoints(:,3);
            seedPoints = seedPoints(:,[2,1]);
            objNum = interp2mex(L,seedPoints(:,1)/params(pi,2),seedPoints(:,2)/params(pi,2));
            objNum = ceil(objNum);
            seedPoints(objNum==0,:) = [];
            objNum(objNum==0) = [];
            
            for drNo = 1:numel(dr)
                idx = rangesearch(NS,seedPoints,drs(drNo));
                TP(ind,ni,pi,drNo) = numel(unique(cat(2,idx{:})));
            end

            TPplusFP(ind,ni,pi) = size(seedPoints,1);

            %Number of centers difference
            numCents = accumarray(objNum,1,[max(objNumbers),1],[],0);
            d = numCents(objNumbers) - correctNumCents(objNumbers);

            d(abs(d)>3) = [];
            d = d + 4;
            dN(ind,ni,pi,:) = accumarray(d,1,[7,1],[],0);
        end
    end
    
    options.Potential_Extent = options.Potential_Extent - 0.5;
    
end
% profile viewer
% save([pth 'resultsInfo_objectScale_ra_sweep_scaleInvarDepth19max_full_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

%
% Analyze and save results ----------------------------------------------
% Results for each individual image
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
dims.PnnRnnF1nn = {'image','trial','sweep','dr'};
dims.PRF1 = {'1','trial','sweep','dr'};
dims.dN = {'image','trial','sweep','dNm3top3'};
dims.TP = {'image','trial','sweep','dr'};
dims.TPplusFPorFN = {'image','trial','sweep'};
dims.dr = dr;
dims.dN = -3:3;
dims.params = params;
dims.paramsDisp = {'ra','scale','r0'};
dims.image = {'LD2P24','LD3P24','LD4P24','LD5P24','LD67P24'};
results.dims = dims;
results.notes = 'History_Size = 5; no speed limit; converge @ 2e-3; no intensity; radLab computer; depth 19 at max DT';
%
% save([pth 'results_objectScale_ra_sweep_scaleInvarDepth19max_full_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'results')
% fprintf('ra = 30, scale = 25 @ 90%%, dapi 67\n')
disp(squeeze(mean(results.F1,2))')
disp(squeeze(mean(sum(results.dN,1),2))'/(numel(names)*484))
disp(squeeze(mean(dN(:,:,:,4),2))'/484)
% figure(1)
% clf
% imshow(I)
% hold on
% plot(seedPoints(:,1),seedPoints(:,2),'.r')
