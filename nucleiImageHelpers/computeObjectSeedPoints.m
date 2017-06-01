function [seedPoints, Info] = computeObjectSeedPoints(BW, M, r0set, useCentroid, options, objNumber, errorCount)
% COMPUTEOBJECTSEEDPOINTS Compute the seed points for a single object
%
% [seedPoints, Info] = computeObjectSeedPoints(BW, M, r0set, useCentroid, options, objNumber)
%
% Input parameters:
%   BW : binary mask of a single object
%   M : base potential modifier, should be 1 outside of object
%   r0set : Nx2 array, set of initial positions
%   useCentroid : logical, if true, then the centroid of the mask is
%       returned.
%   objNumber : Optional. The object number of this object (this is used
%       for creating error messages only)
%
% Output parameters:
%   seedPoints : Lx2 array of computed seed points
%   Info : If options.debug is False, then Info will be empty. If
%       options.debug is True, then Info will be a structure with several
%       fields:
%           r0 : initial points used in simulation
%           r_final : final location of simulated points
%           V : confining potential
%           ComputeInitialPoints : structure of initial points information,
%               see computeInitialPoints()
%           message : integer giving an exit code
%               0, everything is fine
%               1, object is convex or too small (the center of the object
%                   will be the seed point)
%               2, less than 2 initial particles (the center of the object
%                   will be the seed point)
%               3, there was twice an error, object will be skipped

% See also MODELPARTICLEDYNAMICS EXTRACTCLUSTERCENTERS COMPUTEINITIALPOINTS
% CREATE_SCALEINVAR_CONFINING_POTENTIAL COMPUTEGRIDSCALEFACTORS

% James Kapaldo

% M : base potential modifier - should be 1 outside of nuclei.


Info = [];

% Get the error count
if nargin < 6
    objNumber = 1;
end
if nargin < 7
    errorCount = 0;
end

DEBUG = options.Debug;

if useCentroid
    [seedPoints, Info] = computeCentroid(BW,DEBUG,1);%'convex_or_tooSmall');
    return
end % if

try
    PAD_SIZE = options.Potential_Padding_Size;

    % Setup problem ------------------------------------------------------
    % Add in potential modifier
    if ~isempty(M)
        options.Potential_Modifier = @(V) V .* M;
    end

    % Pad initial positions
    r0set = r0set + PAD_SIZE;

    % Compute confining force, initial points, and problem scales.
    [dV, r0, problem_scales, SetupInfo] = setup_problem(BW, [], options, r0set);
    solver_to_data = @(x) (x./problem_scales.grid_to_solver - PAD_SIZE);

    % Model dynamics -----------------------------------------------------
    R = options.Iterations;
    seedPoints_n = cell(R,1);

    if DEBUG
        Info = struct('simulationInfo', cell(1,R), ...
                      'r0', cell(1,R), ...
                      'r_final', cell(1,R), ...
                      'seedPoints_n', cell(1,R), ...
                      'cluster_sizes_n', cell(1,R), ...
                      'iteration_sizes', cell(1,R));
        for n = 1:R
            [r, simInfo] = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
            [seedPoints_n{n}, clstSz] = extractClusterCenters(r, options);

            Info.simulationInfo{n} = simInfo;
            Info.r0{n} = solver_to_data(r0);
            Info.r_final{n} = solver_to_data(r);
            Info.seedPoints_n{n} = solver_to_data(seedPoints_n{n});
            Info.cluster_sizes_n{n} = clstSz;
            Info.iteration_sizes{n} = [size(r,1), numel(clstSz)];
        end
    else
        for n = 1:R
            r = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
            [seedPoints_n{n}, clstSz] = extractClusterCenters(r, options);
        end
    end

    % % Post-processing --------------------------------------------------
    % % Remove any particles outide the object.
    seedPoint_set = cat(1, seedPoints_n{:});
    % Vi = interpolatePotential(V, seedPoint_set, problem_scales);
    % seedPoint_set(Vi>1,:) = [];

    % Cluster the results of the iterations
    if R > 1
        [seedPoints, clstSz] = extractClusterCenters(seedPoint_set, options);
    else
        seedPoints = seedPoint_set;
    end
    toRemove = clstSz < options.Minimum_Cluster_Size;
    seedPoints(toRemove,:) = [];
    seedPoints = solver_to_data(seedPoints);

    if DEBUG
        Info.cluster_sizes = clstSz;
        Info.V = SetupInfo.V;
        Info.ComputeInitialPoints = SetupInfo.ComputeInitialPointsInfo;

        to_combine = {'r0','r_final','seedPoints_n','cluster_sizes_n','iteration_sizes'};
        for fn = 1:numel(to_combine)
            Info.(to_combine{fn}) = cat(1, Info.(to_combine{fn}){:});
        end
        Info.solverTime = mean(cellfun(@(x) x.solverTime, Info.simulationInfo))
        Info.message = 0;%'';
    end

catch ME
    % It can be somtimes that there is an error due to a particle flying out of the image region, or something that rarely (hopefully) happens; therefore, we will try to run the function again. If there is an error a second time, then we issue a warning, save the error, and continue on to the next object.
    if errorCount < 1
        [seedPoints, Info] = computeObjectSeedPoints(BW, M, r0set, useCentroid, options, objNumber, errorCount+1);
    else
        seedPoints = [NaN, NaN];
        if DEBUG
            Info = emptyInfo();
            Info.message = 3;
            Info.error = ME;
        else
            Info.error = ME;
        end
        fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', objNumber, objNumber)
        fprintf(2,'%s\n', getReport(Info.error,'basic'))
        fprintf('\n')
    end
end

end

function [seedPoint, Info] = computeCentroid(BW,DEBUG,reason)
[i,j] = find(BW);
seedPoint = mean([i,j],1);
Info = [];
if DEBUG
    Info = emptyInfo();
    Info.message = reason;
end
end

% function Vi = interpolatePotential(V, r, problem_scales)
%
% r = r ./ problem_scales.grid_to_solver;
% Vi = interp2mex(V, r(:,2), r(:,1))
%
% end


function Info = emptyInfo()
    Info = struct('ComputeInitialPoints', struct(), ...
                  'simulationInfo', {struct()}, ...
                  'r0', NaN(1,2), ...
                  'r_final', NaN(1,2), ...
                  'seedPoints_n', NaN(1,2), ...
                  'cluster_sizes_n', NaN, ...
                  'iteration_sizes', NaN(1,2), ...
                  'cluster_sizes', NaN, ...
                  'V', [], ...
                  'solverTime', NaN, ...
                  'message', NaN);
end
