%% Test restructured code
pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];


% Initialize options ----------------------------------------------------

WIGNER_SEITZ_RADIUS = 5;
CURVATURE_MAX_RADIUS = 35;
CURVATURE_SMOOTHING_SIZE = 2;
INITIAL_SPEED = 0.1;
MINIMUM_HOLE_SIZE = 15;

options = seedPointOptions();
options.Wigner_Seitz_Radius         = 5;
options.Potential_Depth             = -1;
options.Potential_Minimum_Location  = 2;
options.Potential_Extent            = 15;
options.Minimum_Hole_Size           = 15;
options.Use_GPU                     = true;
options.Use_Parallel                = true;
options.Point_Selection_Method      = 'r0set_uniformRandom';
options.Debug                       = true;
options.Initial_Speed               = INITIAL_SPEED;
options.Potential_Scale             = 7;
% options.Maximum_Initial_Potential   = 1/2;
% options.Object_Of_Interest          = 152;


% Set up test parameters ------------------------------------------------
names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10;
N = 5;


% R = [10:5:40, 50:10:100];

% Already computed ===
% R = 10:5:40;
% th = 45:5:85;
% scales = 1;
% [R, Th, scale] = meshgrid(R,th,scales);
% R(1:2:end,:,:) = R(1:2:end,:,:) - 2.5*(1 + (R(1:2:end,:,:)>40));
% ra = R.*sind(Th);
% ps = R.*cosd(Th);
% paramsDone = [ra(:),ps(:),round(ra(:)*0.12),scale(:)];

% new to compute ===
% R = 17:1:34;
% th= 30:2:48;
% scales = 1;
% 
% [R,Th,scale] = meshgrid(R,th,scales);
% % R(1:2:end,:,:) = R(1:2:end,:,:) - 2.5*(1 + (R(1:2:end,:,:)>40));
% R(1:2:end,:) = R(1:2:end,:) + 0.5;

[ra, ps, scale, speed] = ndgrid(13:15, 18:20, 1, 0.01);

ra = ra(:);
ps = ps(:);
scale = scale(:);

% R = sqrt(ra.^2 + ps.^2);
% Th = atand(ra./ps);
% 
% bad = R < 17 | R > 34 | Th < 30 | Th > 48;
% ra(bad) = [];
% ps(bad) = [];
% scale(bad) = [];

% R = R-2.5;
% ra = R.*sind(Th);
% ps = R.*cosd(Th);
% 
% 
% bad = ra<5 | ps < 3 | ps >28 | ra > 25;
% ps(bad) = [];
% ra(bad) = [];
% scale(bad) = [];

r0 = round(ra*0.12);
% params = [ra(:),ps(:),r0(:),scale(:),speed(:)];
params = [14,19,2,1,0.01];
% params(ismembertol(params,paramsDone,'ByRows',0.3),:) = [];
% params = round(params,1);

% figure
% hold on
% try delete(findall(get(gca,'Children'),'Tag','testPoints')), catch, end
% line(params(:,2),params(:,1),params(:,4),'Marker','.','Color','r','linestyle','none','tag','testPoints')
% daspect([1 1 1])

% [ra,ps,scale] = meshgrid(10:2.5:25, 7, [0.5,1,2,3]);
% r0 = round(ra*0.12);
% params = [ra(:),ps(:),r0(:),scale(:)];


% scale = 1;

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
    TPplusFN(ind,:) = trueNumberOfNuclei;

    for pi = 1:size(params,1)
        if ~mod(pi,2)
            fprintf('  %s >> param %d/%d...\n', datestr(now,31),pi,size(params,1))
        end
        
        if (strcmp(names{ind},'67') || strcmp(names{ind},'5')) && (params(pi,4) > 4)
            options.Use_Parallel = false;
        else
            options.Use_Parallel = true;
        end
        
        options.Potential_Minimum_Location = params(pi,3);
        options.Potential_Extent = params(pi,1);
        options.Potential_Scale = params(pi,2);
        options.Initial_Speed = params(pi,5);
        
        options.Wigner_Seitz_Radius = round( params(pi,4) * WIGNER_SEITZ_RADIUS);
%         options.Curvature_Max_Radius = ceil( params(pi,4) * CURVATURE_MAX_RADIUS);
%         options.Curvature_Smoothing_Size = max(1, ceil( params(pi,4) * CURVATURE_SMOOTHING_SIZE));
        
        options.Minimum_Hole_Size = round( params(pi,4) * MINIMUM_HOLE_SIZE);

        BWs = imresize(BW,params(pi,4),'bilinear','Antialiasing',false)>0.5;
        NS = createns(truthDat(:,[3,2])*params(pi,4));
        drs = dr*params(pi,4);
        
        for ni = 1:N
            if ~mod(ni,10)
                fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
            end
    %         profile on
            start = tic;
            [seedPoints, runInfo] = computeNucleiCenters_distTransform(I,BWs,options);
            totalTime = toc(start);
    %         profile off
            

            % If we are scaling the objects, then the number of objects can
            % change due to objects joining or seperating. We need all of
            % the object numbers to correspond to the original objects;
            % thus, for each seed point we compute what object is is in.
            if params(pi,4) ~= 1
                seedPoints(:,3) = interp2mex(L,seedPoints(:,2)/params(pi,4),seedPoints(:,1)/params(pi,4));
                seedPoints(:,3) = ceil(seedPoints(:,3));
                seedPoints(seedPoints(:,3)==0,:) = [];
            end

            if options.Debug
                Info(ind,ni,pi).centers = seedPoints;
                Info(ind,ni,pi).solverTime = cellfun(@(x) x.solverTime, runInfo);
                Info(ind,ni,pi).totalComputationTime = totalTime;
                Info(ind,ni,pi).N = cellfun(@(x) size(x.r0,1), runInfo);
                Info(ind,ni,pi).message = cellfun(@(x) x.message, runInfo);
            end
            
            clear runInfo
            
            objNum = seedPoints(:,3);
            seedPoints = seedPoints(:,[2,1]);
            
            % Compute TP+FP: this is the number of seed points calculated
            % for each object
            TPplusFP(ind, ni, pi) = size(seedPoints,1);
            
            % Compute TP
            [nnIdx,nnD] = knnsearch(NS,seedPoints);
            
            toRemove = nnD >= drs(end);
            nnD(toRemove) = [];
            nnIdx(toRemove) = [];
            objNum(toRemove) = [];
            
            for drNo = 1:numel(dr)
                inRange = nnD < drs(drNo);
                TP(ind, ni, pi, drNo) = numel(unique(nnIdx(inRange)));
            end

            %Number of centers difference
            numCents = accumarray(objNum,1,[max(objNumbers),1],[],0);
            d = numCents(objNumbers) - correctNumCents(objNumbers);

            d(abs(d)>3) = [];
            d = d + 4;
            dN(ind,ni,pi,:) = accumarray(d,1,[7,1],[],0);
        end
    end
end
% profile viewer
save([pth 'rI_potDatMax_ra_speed1_test_dense2_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

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
dims.paramsDisp = {'ra','depth','r0','scale','speed'};
dims.image = {'LD2P24','LD3P24','LD4P24','LD5P24','LD67P24'};
results.dims = dims;
results.notes = 'History_Size = 5; no speed limit; converge @ 2e-3; no intensity; potential depth set at max bwdist; scaling potential and using the scale factor for the pixel size in the gradient calculation. Scale all particle positions (not just the distance between the particles), and unscale when evaluating the potential gradient.';

save([pth 'r_potDatMax_ra_speed1_test_dense2_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'results')
% fprintf('ra = 30, scale = 25 @ 90%%, dapi 67\n')
disp(squeeze(mean(results.F1,2)))
disp(squeeze(mean(sum(results.dN,1),2))/(numel(names)*484))

% figure(1)
% clf
% imshow(I)
% hold on
% plot(seedPoints(:,1),seedPoints(:,2),'.r')
