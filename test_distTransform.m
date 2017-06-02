%% Test restructured code
pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];


% Initialize options ----------------------------------------------------

options = seedPointOptions();
options.Wigner_Seitz_Radius         = 5;
options.Potential_Depth             = -1;
options.Potential_Minimum_Location  = 2;
options.Potential_Extent            = 12;
options.Minimum_Hole_Size           = 15;
options.Point_Selection_Method      = 'r0set_uniformRandom';
options.Use_GPU                     = false;
options.Use_Parallel                = true;
options.Debug                       = false;

options.Potential_Scale             = 17;
% options.Mass_Charge_Multiplier      = 0.1;%19 / options.Potential_Scale;

% options.Maximum_Initial_Potential   = 1/3;
% options.Object_Of_Interest          = 152;

% options.Particle_Damping_Rate = 5e-2;

% Set up test parameters ------------------------------------------------
names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10;
N = 1;

% Initialize variables for storing resuts -------------------------------
TP = zeros(n,N,numel(dr));
TPplusFP = zeros(n,N);
TPplusFN = zeros(n,N);
dN = zeros(n,N,7);

% Info = struct('centers',[],'solverTime',[],'totalComputationTime',NaN,'N',NaN,'message',[]);
% Info(n,N).totalComputationTime = NaN;
% Info = cell(484,1);

% Compute results -------------------------------------------------------

for ind = n:-1:1
    %     options.Potential_Extent            = 15 - 0.5*ind;
    fprintf('%s >> Image = LD%sP24 (%d/%d)...\n', datestr(now,31),names{ind},ind,numel(names))
    
    I = imread(im_pth(names{ind}));
    BW = imread(bw_pth(names{ind})) > 0;
    
    results = load(results_pth(names{ind}));
    result_fields = fieldnames(results);
    tmp = results.(result_fields{1});
    
    NS = createns(tmp(:,[3,2]));
    correctNumCents = accumarray(tmp(:,1),1,[],[],0);
    objNumbers = unique(tmp(:,1));
    trueNumberOfNuclei = size(tmp,1);
    TPplusFN(ind,:) = trueNumberOfNuclei;
    
    for ni = 1:N
        if ~mod(ni,10)
            fprintf('  %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
        end
        %         profile on
        start = tic;
        [seedPoints, Info] = computeNucleiCenters(I,BW,options);
        totalTime = toc(start);
        %         profile off
        %         Info(ind,ni).centers = seedPoints;
        %         Info(ind,ni).solverTime = cellfun(@(x) x.solverTime, runInfo);
        %         Info(ind,ni).totalComputationTime = totalTime;
        %         Info(ind,ni).N = cellfun(@(x) size(x.r0,1), runInfo);
        %         Info(ind,ni).message = cellfun(@(x) x.message, runInfo);
        
        objNum = seedPoints(:,3);
        seedPoints = seedPoints(:,[2,1]);
        
        % Compute TP+FP: this is the number of seed points calculated
        % for each object
        TPplusFP(ind, ni) = size(seedPoints,1);
        
        %Number of centers difference
        numCents = accumarray(objNum,1,[max(objNumbers),1],[],0);
        d = numCents(objNumbers) - correctNumCents(objNumbers);
        
        d(abs(d)>3) = [];
        d = d + 4;
        dN(ind,ni,:) = accumarray(d,1,[7,1],[],0);
        
        % Compute TP
        [nnIdx,nnD] = knnsearch(NS,seedPoints);
        
        toRemove = nnD >= dr(end);
        nnD(toRemove) = [];
        nnIdx(toRemove) = [];
        objNum(toRemove) = [];
        
        for drNo = 1:numel(dr)
            inRange = nnD < dr(drNo);
            TP(ind, ni, drNo) = numel(unique(nnIdx(inRange)));
        end

    end
end
% profile viewer
% save([pth 'resultsInfo_distTransform_standard_N20_v20170120-reStructured_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

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
dims.PnnRnnF1nn = {'image','trial','dr'};
dims.PRF1 = {'1','trial','dr'};
dims.dN = {'image','trial','dNm3top3'};
dims.TP = {'image','trial','dr'};
dims.TPplusFPorFN = {'image','trial'};
dims.dr = dr;
dims.dN = -3:3;
dims.image = {'LD2P24','LD3P24','LD4P24','LD5P24','LD67P24'};
results.dims = dims;

% save([pth 'results_distTransform_standard_N20_v20170120-reStructured' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'results')
% fprintf('ra = 12, scale = 12\n')
disp(squeeze(mean(results.F1,2))')
disp(squeeze(mean(sum(results.dN,1),2))'/(numel(names)*484))
disp(squeeze(mean(results.dN(:,:,3:5),2)/484))
% figure(1)
% clf
% imshow(I)
% hold on
% plot(seedPoints(:,1),seedPoints(:,2),'.r')
