%% Test restructured code
pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];


% Initialize options ----------------------------------------------------

options = seedPointOptions();
options.Max_Distance_Transform      = 18; % This makes it scale invariant.
options.Use_GPU                     = true;
options.Use_Parallel                = true;
options.Debug                       = true;

WIGNER_SEITZ_RADIUS = 5;
MINIMUM_HOLE_SIZE = 15;
POTENTIAL_DEPTH = -1;

% Set up test parameters ------------------------------------------------
names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10; 
N = 1; % number of iterations

% Generate grid of scan points.
[scales,ra] = meshgrid(0.7:0.3:4.9, 6:4:34);
ra(:,1:2:end) = ra(:,1:2:end) + 2;
ra = ra(:);
scales = scales(:);

bad = ra > 35;
ra(bad) = [];
scales(bad) = [];

r0 = round(ra*0.12);
params = [ra,scales,r0];

% figure
% line(params(:,2),params(:,1),ones(size(params,1),1),'Marker','.','Color','r','linestyle','none','tag','testPoints')

% Initialize variables for storing resuts -------------------------------
TP = zeros(n,N,size(params,1),numel(dr));
TPplusFP = zeros(n,N,size(params,1));
TPplusFN = zeros(n,N,size(params,1));
dN = zeros(n,N,size(params,1),7);

Info = struct('centers',[],'solverTime',[],'totalComputationTime',NaN,'N',NaN,'message',[]);
Info(n,N,size(params,1)).totalComputationTime = NaN;

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
    TPplusFN(ind,:) = size(tmp,1); % true Number Of Nuclei

    for pp = 1:size(params,1)
        if ~mod(pp,2)
            fprintf('  %s >> param %d/%d...\n', datestr(now,31),pp,size(params,1))
        end
        
        if (strcmp(names{ind},'67') || strcmp(names{ind},'5')) && (params(pp,2) > 4)
            options.Use_Parallel = false;
        else
            options.Use_Parallel = true;
        end
        
        options.Potential_Parameters = [POTENTIAL_DEPTH, params(pp,[3,1])];
        options.ScaleInvarient_Potential_Extent = params(pp,1);
        options.Wigner_Seitz_Radius = round( params(pp,2) * WIGNER_SEITZ_RADIUS);
        options.Minimum_Hole_Size = round( params(pp,2) * MINIMUM_HOLE_SIZE);

        BWs = imresize(BW,params(pp,2),'bilinear','Antialiasing',false)>0.5;
        NS = createns(truthDat(:,[3,2])*params(pp,2));
        drs = dr*params(pp,2);
        
        for ni = 1:N
            if ~mod(ni,10)
                fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
            end

            start = tic;
            [seedPoints, runInfo] = computeNucleiCenters_distTransform(I,BWs,options);
            totalTime = toc(start);

            % If scaling the objects, then the number of objects can change
            % due to objects joining or seperating. We need all of the
            % object numbers to correspond to the original objects; thus,
            % for each seed point we compute what object it is in.
            if params(pp,2) ~= 1
                seedPoints(:,3) = interp2mex(L,seedPoints(:,2)/params(pp,2),seedPoints(:,1)/params(pp,2));
                seedPoints(:,3) = ceil(seedPoints(:,3));
                seedPoints(seedPoints(:,3)==0,:) = [];
            end

            if options.Debug
                Info(ind,ni,pp).centers = seedPoints;
                Info(ind,ni,pp).solverTime = cellfun(@(x) x.solverTime, runInfo);
                Info(ind,ni,pp).totalComputationTime = totalTime;
                Info(ind,ni,pp).N = cellfun(@(x) size(x.r0,1), runInfo);
                Info(ind,ni,pp).message = cellfun(@(x) x.message, runInfo);
            end
            
            clear runInfo
            
            objNum = seedPoints(:,3);
            seedPoints = seedPoints(:,[2,1]);
            
            % Compute TP+FP: this is the number of seed points calculated
            % for each object
            TPplusFP(ind, ni, pp) = size(seedPoints,1);
            
            %Number of centers difference
            numCents = accumarray(objNum,1,[max(objNumbers),1],[],0);
            d = numCents(objNumbers) - correctNumCents(objNumbers);

            d(abs(d)>3) = [];
            d = d + 4;
            dN(ind,ni,pp,:) = accumarray(d,1,[7,1],[],0);
            
            % Compute TP
            [nnIdx,nnD] = knnsearch(NS,seedPoints);
            
            toRemove = nnD >= drs(end);
            nnD(toRemove) = [];
            nnIdx(toRemove) = [];
            objNum(toRemove) = [];
            
            for drNo = 1:numel(dr)
                inRange = nnD < drs(drNo);
                TP(ind, ni, pp, drNo) = numel(unique(nnIdx(inRange)));
            end
        end
    end
end

save([pth 'resultsInfo_scale_ra_scan_scaleInvar_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

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
results.notes = 'History_Size = 5; no speed limit; converge @ 2e-3; commit e3ee7460cc13dd1c21b54396c45b0746fa33e4d4 (Tue Feb 7 21:24:22 2017 -0500)';
results.options = options;

save([pth 'results_scale_ra_scan_scaleInvar_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'results')