function [truth, nuclei_per_object] = load_truth(trial_str)
% LOAD_TRUTH Load truth data about the nuclei center locations.
%
% [truth, nuclei_per_object] = load_truth(trial_str)
%
% Input parameters:
% trial_str : A string that can have the value of {'2','3','4','5','67'}.
%
% Output parameters: 
% truth : An Nx3 array where each row gives the true location of a single
%   nuclei. The first column gives the index of the object containing the
%   nuclei. The second and third columns give the y,x location of the
%   nuclei.
% nuclei_per_object: Mx1 array, where M is the number of objects, giving
%   the number of nuclei per object.

% James Kapaldo

    if ~contains(trial_str, {'2','3','4','5','67'})
        error('load_truth:unknownTrial','Unknown trial. Allowed trials are ''2'', ''3'', ''4'', ''5'', or ''67''.')
    end
    
    truth = load(['exampleData\markedCenters_' trial_str '']);
    truth_fields = fieldnames(truth);
    truth = truth.(truth_fields{1});
    
    nuclei_per_object = accumarray(truth(:,1),1,[],[],0);
end


%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
