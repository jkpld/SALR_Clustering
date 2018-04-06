function [F1, dN] = evaluatePerformance(nucleiCenters, dr)
% EVALUATEPERFORMANCE Compute the F1 and fractional distribution of the
% difference between the true number of nuclei in each clump and the
% computed number of nuclei in each clump.
%
% [F1, dN] = evaluatePerformance(nucleiCenters, dr)
%
% Input parameters:
% nucleiCenters : 1x5 cell array containing the cmoputed nuclei centers for
%   each of the 5 images, as returned by computeNucleiCenters. They are
%   assumed to be in the order of '2','3','4','5','67'.
% dr : An array of distances at which to evaluate the F1 score.

% James Kapaldo

trials = {'2','3','4','5','67'};

TP = zeros(1,numel(dr)); % Number of true positives
TPplusFP = 0; % Total number of computed nuclei
TPplusFN = 0; % Total number of true nuclei
dN = zeros(1,7); % The fractional distribution of the difference between the true number of nuclei and the computed number of nuclei.

for t = 1:numel(trials)
    % The length of nucleiCenters gives the number of nuclei found.
    TPplusFP = TPplusFP + size(nucleiCenters{t},1);
    
    % Load the truth data. truth contains the location of the true nuclei
    % and nuclei_per_object contains the number of nuclei in each clump
    [truth, nuclei_per_object] = load_truth(trials{t});
    
    % The number of true nuclei
    TPplusFN = TPplusFN + size(truth,1);
    
    % Compute the difference between the true number of nuclei in each
    % object and the computed number of nuclei in each object. The index of
    % the object to which a nuclei center belongs is in the third column of
    % nucleiCenters.
    n = accumarray(nucleiCenters{t}(:,3),1,[],[],0); % number of nuclei in each clump
    d = n - nuclei_per_object; % difference
    
    % Compute the fraction of nuclei clumps with d=-3 to d=+3.
    d(abs(d)>3) = [];
    dN = dN + accumarray(d+4,1,[7,1],[],0)/numel(d)/numel(trials);
    
    % Create a nearest neighbor searcher
    NS = createns(truth(:,[3,2]));
    
    % Compute TP
    [nnIdx,nnD] = knnsearch(NS,nucleiCenters{t}(:,[2,1]));

    for drNo = 1:numel(dr)
        inRange = nnD < dr(drNo);
        TP(drNo) = TP(drNo) + numel(unique(nnIdx(inRange)));
    end
end

P = TP./TPplusFP; % precision
R = TP./TPplusFN; % recall
F1 = 2*P.*R./(P+R); % F1 score

end

%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
