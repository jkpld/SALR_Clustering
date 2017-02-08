%% Test restructured code
pth = 'K:\Google_Drive\MATLAB\seed_point_detection\';

im_pth      = @(n) [pth 'exampleImages\testImage_image_LD' n 'P24.tif'];
bw_pth      = @(n) [pth 'exampleImages\testImage_mask_LD' n 'P24.tif'];
results_pth = @(n) [pth 'exampleImages\markedCenters_LD' n 'P24'];


% Initialize options ----------------------------------------------------

WIGNER_SEITZ_RADIUS = 5;
CURVATURE_MAX_RADIUS = 35;
CURVATURE_SMOOTHING_SIZE = 2;
MINIMUM_HOLE_SIZE = 15;

options = seedPointOptions();
options.Wigner_Seitz_Radius         = WIGNER_SEITZ_RADIUS;
options.Potential_Depth             = -1;
options.Potential_Minimum_Location  = 2;
options.Potential_Extent            = 15;
options.Minimum_Hole_Size           = MINIMUM_HOLE_SIZE;
options.Point_Selection_Method      = 'r0set_uniformRandom';
options.Use_GPU                     = true;
options.Use_Parallel                = true;
options.Debug                       = false;

options.Potential_Scale             = 19;
% options.Object_Of_Interest = 125;
% Set up test parameters ------------------------------------------------
names = {'2','3','4','5','67'};
n = numel(names);
dr = 3:10;
N = 1;

ra = 12:20;
ps = 18;
scale = 1;

[ra,ps,scale] = meshgrid(ra,ps,scale);

r0 = round(ra*0.12);
params = [ra(:),ps(:),r0(:),scale(:)];


% Initialize variables for storing resuts -------------------------------
% Want to store the results for each object. There are 484 objects per
% image.


TP          = zeros(n, 484, N, size(params,1), numel(dr));
TPplusFP    = zeros(n, 484, N, size(params,1));
TPplusFN    = zeros(n, 484);

TP2          = zeros(n, N, size(params,1), numel(dr));
TPplusFP2    = zeros(n, N, size(params,1));
TPplusFN2    = zeros(n, 1);

% Info = struct('centers',[],'solverTime',[],'totalComputationTime',NaN,'N',NaN,'message',[]);
% Info(n,N,size(params,1)).totalComputationTime = NaN;

% Compute results -------------------------------------------------------
dr = sort(dr,'ascend');
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
    TPplusFN(ind,:) = correctNumCents;
    
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
        options.Potential_Scale = params(pi,2);
        
        options.Wigner_Seitz_Radius         = round( params(pi,4) * WIGNER_SEITZ_RADIUS);
        options.Curvature_Max_Radius        = ceil( params(pi,4) * CURVATURE_MAX_RADIUS);
        options.Curvature_Smoothing_Size    = max(1, ceil( params(pi,4) * CURVATURE_SMOOTHING_SIZE));
        
        options.Minimum_Hole_Size = round( params(pi,4) * MINIMUM_HOLE_SIZE);

        BWs = imresize(BW,params(pi,4),'bilinear','Antialiasing',false)>0.5;
        NS = createns(truthDat(:,[3,2])*params(pi,4));
        drs = dr*params(pi,4);
        
        for ni = 1:N
            if ~mod(ni,10)
                fprintf('    %s >> iteration %d/%d...\n', datestr(now,31),ni,N)
            end
    
            start = tic;
            [seedPoints, runInfo] = computeNucleiCenters_distTransform(I,BWs,options);
            totalTime = toc(start);

            % If we are scaling the objects, then the number of objects can
            % change due to objects joining or seperating. We need all of
            % the object numbers to correspond to the original objects;
            % thus, for each seed point we compute what object is is in.
            if params(pi,4) ~= 1
                seedPoints(:,3) = interp2mex(L,seedPoints(:,2)/params(pi,4),seedPoints(:,1)/params(pi,4));
                seedPoints(:,3) = ceil(seedPoints(:,3));
                seedPoints(seedPoints(:,3)==0,:) = [];
            end
            
            % Save some info results
%             Info(ind,ni,pi).centers = seedPoints;
%             Info(ind,ni,pi).solverTime = cellfun(@(x) x.solverTime, runInfo);
%             Info(ind,ni,pi).totalComputationTime = totalTime;
%             Info(ind,ni,pi).N = cellfun(@(x) size(x.r0,1), runInfo);
%             Info(ind,ni,pi).message = cellfun(@(x) x.message, runInfo);

            clear runInfo
            
            objNum = seedPoints(:,3);
            seedPoints = seedPoints(:,[2,1]);
            
            % Compute TP+FP: this is the number of seed points calculated
            % for each object
            TPplusFP(ind, :, ni, pi) = accumarray(objNum, 1, [max(objNumbers),1], [], 0);
            
            % Compute TP
            [nnIdx,nnD] = knnsearch(NS,seedPoints);
            
            toRemove = nnD >= drs(end);
            nnD(toRemove) = [];
            nnIdx(toRemove) = [];
            objNum(toRemove) = [];
            
            for drNo = 1:numel(dr)
                inRange = nnD < drs(drNo);
                TP(ind, :, ni, pi, drNo) = accumarray(objNum(inRange), nnIdx(inRange), [max(objNumbers), 1], @(x) numel(unique(x)), 0);
            end
        end
    end
end
% profile viewer
% save([pth 'resultsInfo_objectScale_ra_sweep_scaleInvarDepth19max_full_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'Info','options')

%
% Analyze and save results ----------------------------------------------
% Results for each individual image
%%
% TP2 = cat(3,TP,TP);R
[TPb, TPplusFPb, TPplusFNb, idx_best_performance,score] = maximizeObjectPerformance(TP, TPplusFP, TPplusFN, 4);

%% create feature and output matrices

params_of_best_performance = params(idx_best_performance(:), 1); % output
features = [round(objMaxDT(:)), objArea(:), objMeanDT(:), objSumPosCurv(:)];


figure
scatter3(features(:,1), features(:,2), features(:,3), 3*reshape(ones(5,484).*(1:5)',5*484,1), params_of_best_performance,'filled');%'marker','.','linestyle','none','color','k')
xlabel('max dt')
ylabel('area/max^2')
zlabel('mean')
drawnow;
goDark(gcf)

%%
figure
line(features(:,3), params_of_best_performance,'marker','.','linestyle','none')

%%
F1 = @(tp, tpfn, tpfp) 2*tp ./ (tpfn + tpfp);

TP3 = TP(:,:,:,:,1);

TP3 = reshape(TP3,[5*484,9]);
tppfp = reshape(TPplusFP, [5*484, 9]);
tppfn = TPplusFN(:) * ones(1,9);

x = (18 ./ objMaxDT(:)) .* ones(1, size(params,1));
y = params(:,1)' .* ones(5*484, 1);

dtp = decimateData(x,y,TP3,'binSize',[0.05,1],'reductionMethod',@sum,'generatePlots',false);
dtppfp = decimateData(x,y,tppfp,'binSize',[0.05,1],'reductionMethod',@sum,'generatePlots',false);
dtppfn = decimateData(x,y,tppfn,'binSize',[0.05,1],'reductionMethod',@sum,'generatePlots',false);

f1 = F1(dtp.Z, dtppfp.Z, dtppfn.Z);

figure
surface(dtp.X, dtp.Y, f1)
axis tight
xlabel('mean dt * (18/max dt)')
ylabel('ra')

%%

x = objMeanDT(:) .* ones(1, size(params,1));
y = params(:,1)' .* ones(5*484, 1);
z = reshape(score, [5*484,9]);
toRemove = all(z==2,2);
% x(toRemove) = [];
% y(toRemove) = [];
% z(toRemove) = [];
%
decimateData(x,y,z,'binSize',[0.5,1],'reductionMethod',@sum,'generatePlots',true);


%%

computeF1 = @(tp, tpfn, tpfp) 2*tp ./ (tpfn + tpfp);
disp here
F1img = squeeze(computeF1(sum(TP,2), sum(TPplusFP,2), sum(TPplusFN,2)))
F1 = squeeze(computeF1(sum(sum(TP,2),1), sum(sum(TPplusFP,2),1), sum(sum(TPplusFN,2),1)))'

dN = squeeze(TPplusFP - TPplusFN);
dN = squeeze(sum(sum(dN==0,1),2)/(5*484))
%%


computeF1 = @(p,r) 2*p.*r./(p+r);

Pobj = TP./TPplusFP;
Robj = TP./TPplusFN;
F1obj = computeF1(Pobj,Robj);




[~,ind] = max(F1obj(:,:,:,:,1),[],4);

tp3 = TP(:,:,:,:,1);

Pimg = sum(tp3(ind),2)./sum(TPplusFP(ind),2);
Rimg = sum(tp3(ind),2)./sum(TPplusFN,2);
F1img = computeF1(Pimg, Rimg);

P = sum(sum(tp3(ind),2),1)./sum(sum(TPplusFP(ind),2),1);
R = sum(sum(tp3(ind),2),1)./sum(sum(TPplusFN,2),1);
F1 = computeF1(P, R);

%%

Pobj = TP./TPplusFP;
Robj = TP./TPplusFN;
F1obj = 2*Pobj.*Robj./(Pobj+Robj);



Pnn2 = TP2./TPplusFP2;
Rnn2 = TP2./TPplusFN2;
F1nn2 = 2*Pnn2.*Rnn2./(Pnn2+Rnn2);

dN = TPplusFP - TPplusFN;

% Results for all images
% P = sum(TP,1)./sum(TPplusFP,1);
% R = sum(TP,1)./sum(TPplusFN,1);
% F1 = 2*P.*R./(P+R);
% 
% results = [];
% results.Pnn = Pnn;
% results.Rnn = Rnn;
% results.F1nn = F1nn;
% results.P = P;
% results.R = R;
% results.F1 = F1;
% results.dN = dN;
% results.TP = TP;
% results.TPplusFP = TPplusFP;
% results.TPplusFN = TPplusFN;
% 
% dims = [];
% dims.PnnRnnF1nn = {'image','trial','sweep','dr'};
% dims.PRF1 = {'1','trial','sweep','dr'};
% dims.dN = {'image','trial','sweep','dNm3top3'};
% dims.TP = {'image','trial','sweep','dr'};
% dims.TPplusFPorFN = {'image','trial','sweep'};
% dims.dr = dr;
% dims.dN = -3:3;
% dims.params = params;
% dims.paramsDisp = {'ra','scale','r0'};
% dims.image = {'LD2P24','LD3P24','LD4P24','LD5P24','LD67P24'};
% results.dims = dims;
% results.notes = 'History_Size = 5; no speed limit; converge @ 2e-3; no intensity; radLab computer; depth 19 at max DT';
%
% save([pth 'results_objectScale_ra_sweep_scaleInvarDepth19max_full_' datestr(now,'yyyymmddTHHMMSS') '.mat'], 'results')
% fprintf('ra = 30, scale = 25 @ 90%%, dapi 67\n')
% disp(squeeze(mean(results.F1,2))')
% disp(squeeze(mean(sum(results.dN,1),2))'/(numel(names)*484))
% disp(squeeze(dN(:,:,:,4))'/484)
% figure(1)
% clf
% imshow(I)
% hold on
% plot(seedPoints(:,1),seedPoints(:,2),'.r')
