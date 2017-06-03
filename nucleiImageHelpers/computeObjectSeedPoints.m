% function [seedPoints, Info] = computeObjectSeedPoints(binned_data, M, r0set, useCentroid, options, objNumber, errorCount)
function [seedPoints, Info] = computeObjectSeedPoints(binned_data, options, varargin)
% COMPUTEOBJECTSEEDPOINTS Compute the seed points for a single object
%
% [seedPoints, Info] = computeObjectSeedPoints(binned_data, M, r0set, useCentroid, options, objNumber)
%
% Input parameters:
%   binned_data : binned data of a single object
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
%           ComputeInitialPointsInfo : structure of initial points
%               information, see computeInitialPoints()
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

[D, data_limits, r0set, modifier, useCentroid, objNumber, errorCount, Use_Parallel, verbose] = ...
    parse_inputs(binned_data, options, varargin{:});

Info = [];
DEBUG = options.Debug;

if useCentroid
    [seedPoints, Info] = computeCentroid(binned_data,DEBUG,1);%'convex_or_tooSmall');
    return
end

try
    % Setup problem ------------------------------------------------------
    % Add in potential modifier
    if ~isempty(modifier)
        options.Potential_Modifier = @(V) V .* modifier;
    end

    % Compute confining force, initial points, and problem scales.
    [dV, r0, problem_scales, V, SetupInfo] = setup_problem(binned_data, data_limits, options, r0set);
    solver_to_data = @(x) problem_scales.grid_to_data(x./problem_scales.grid_to_solver);

    % Model dynamics -----------------------------------------------------
    % If there was only one initial particle, then we do not simulate it,
    % we will just return the centroid of the object.
    %
    % Note: Since we can model the object several times, if all of the
    % iterations have less than 4 particles and at least 1 of them only has
    % one particle, then we will not simulation it.
    num_r0 = cellfun(@(x) size(x,1), r0);
    if any(num_r0<2) && all(num_r0<4)
        [seedPoints,Info] = computeCentroid(binned_data,DEBUG,2);%'lessThan_2_initial_particles');
        if DEBUG
            Info.ComputeInitialPoints = SetupInfo.ComputeInitialPointsInfo;
        end
        return;
    end

    R = options.Iterations;
    if Use_Parallel && R == 1
        Use_Parallel = false;
    end

    seedPoints_n = cell(R,1);
    cluster_sizes_n = cell(R,1);

    if DEBUG
        r_final_n = cell(R,1);
        iteration_sizes = zeros(R,2);
        simulationInfo = cell(1,R);

        if Use_Parallel
            parfor n = 1:R
                [r_final_n{n}, simulationInfo{n}] = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r_final_n{n}, options);
                iteration_sizes(n,:) = [size(r_final_n{n},1), numel(cluster_sizes_n{n})];
            end
        else
            for n = 1:R
                [r_final_n{n}, simulationInfo{n}] = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r_final_n{n}, options);
                iteration_sizes(n,:) = [size(r_final_n{n},1), numel(cluster_sizes_n{n})];
            end
        end

        Info.simulationInfo = simulationInfo;
        Info.r0 = cellfun(@(x) solver_to_data(x), r0, 'UniformOutput', false);
        Info.r_final = cellfun(@(x) solver_to_data(x), r_final_n, 'UniformOutput', false);
        Info.seedPoints_n = cellfun(@(x) solver_to_data(x), seedPoints_n, 'UniformOutput', false);
        Info.cluster_sizes_n = cluster_sizes_n;
        Info.iteration_sizes = iteration_sizes;

    else
        if Use_Parallel
            parfor n = 1:R
                r = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r, options);
            end
        else
            for n = 1:R
                r = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r, options);
            end
        end
    end

    % % Post-processing --------------------------------------------------
    % Remove any particles outide the object.
    seedPoint_set = cat(1, seedPoints_n{:});
    Vi = interpolatePotential(V, seedPoint_set ./ problem_scales.grid_to_solver);
    seedPoint_set(Vi>1,:) = [];

    % Cluster the results of the iterations
    if R > 1
        [seedPoints, clstSz] = extractClusterCenters(seedPoint_set, options);
    else
        seedPoints = seedPoint_set;
        clstSz = cluster_sizes_n{1};
    end
    toRemove = clstSz < options.Minimum_Cluster_Size;
    seedPoints(toRemove,:) = [];
    seedPoints = solver_to_data(seedPoints);

    if DEBUG
        Info.cluster_sizes = clstSz;
        Info.V = V;
        Info.ComputeInitialPoints = SetupInfo.ComputeInitialPointsInfo;

        to_combine = {'r0','r_final','seedPoints_n','cluster_sizes_n'};
        for fn = 1:numel(to_combine)
            Info.(to_combine{fn}) = cat(1, Info.(to_combine{fn}){:});
        end
        Info.solverTime = mean(cellfun(@(x) x.solverTime, Info.simulationInfo),'omitnan');
        Info.message = 0;%'';
    end

catch ME
    % It can be somtimes that there is an error in modelParticleDynamics
    % due to a particle flying out of the image region, or something that
    % rarely (hopefully) happens; therefore, we will try to run the
    % function again. If there is an error a second time, then we issue a
    % warning, save the error, and continue on to the next object.
    if any(strcmp('modelParticleDynamics',{ME.stack.name}))
        if errorCount < 1
            inputs.data_limits = data_limits;
            inputs.modifier = M;
            inputs.r0set = r0set;
            inputs.useCentroid = useCentroid;
            inputs.objNumber = objNumber;
            inputs.verbose = verbose;
            inputs.errorCount = errorCount + 1;
            [seedPoints, Info] = computeObjectSeedPoints(binned_data, options, inputs);
        else
            seedPoints = NaN(1,D);
            if DEBUG
                Info = emptyInfo(D);
                Info.message = 3;
                Info.error = ME;
            else
                Info.error = ME;
            end
            fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', objNumber, objNumber)
            fprintf(2,'%s\n', getReport(Info.error,'basic'))
            fprintf('\n')
        end
    else
        rethrow(ME)
    end
end

end

function [seedPoint, Info] = computeCentroid(binned_data,DEBUG,reason)

% Note, binned_data has not been padded and so the centroid is already in data
% units.
sz = size(binned_data);
D = numel(sz);
[BWpts{1:D}] = ind2sub(sz,find(binned_data));
BWpts = cat(2,BWpts{:});
seedPoint = mean(BWpts,1);

Info = [];
if DEBUG
    Info = emptyInfo(D);
    Info.message = reason;
end
end

function Vi = interpolatePotential(V, r)

D = size(r,2);

if D == 2
    Vi = interp2mex(V, r(:,2), r(:,1));
else
    Vi = griddedInterpolant(V);
    Vi = Vi(r);
end

end


function Info = emptyInfo(D)
    Info = struct('ComputeInitialPoints', struct(), ...
                  'simulationInfo', {struct('solverTime',NaN,'ode_solution',[],'SystemInputs',[])}, ...
                  'r0', NaN(1,D), ...
                  'r_final', NaN(1,D), ...
                  'seedPoints_n', NaN(1,D), ...
                  'cluster_sizes_n', NaN, ...
                  'iteration_sizes', NaN(1,D), ...
                  'cluster_sizes', NaN, ...
                  'V', [], ...
                  'solverTime', NaN, ...
                  'message', NaN);
end


function [D, data_limits, r0set, modifier, useCentroid, objNumber, errorCount, Use_Parallel, verbose] = parse_inputs(binned_data, options, varargin)

    sz = size(binned_data); % Size
    D = numel(sz); % Dimension
    PAD_SIZE = options.Potential_Padding_Size;

    % Get input values
    p = inputParser;
    p.FunctionName = 'copmuteObjectSeedPoints';
    
    function valid = validate_r0set(t)
        if ~isempty(t)
            validateattributes(t, {'double'}, {'size',[NaN, D],'finite','real'},'varName','r0set')
        end
        valid = true;
    end

    function valid = validate_datalimits(t)
        if ~isempty(t)
            validateattributes(t, {'numeric'}, {'size',[2, D],'finite','real'},'varName','data_limits')
        end
        valid = true;
    end

    function valid = validate_modifier(t)
        if ~isempty(t)
            validateattributes(t, {'double'}, {'size',sz + 2*PAD_SIZE, 'finite', 'real'},'varName','modifier')
        end
        valid = true;
    end

    addParameter(p,'data_limits',[], @validate_datalimits)
    addParameter(p,'r0set',[], @validate_r0set)
    addParameter(p,'modifier',[], @validate_modifier)
    addParameter(p,'useCentroid',false, @(t) t==0 || t==1)
    addParameter(p,'objNumber',1, @(t) validateattributes(t, {'numeric'}, {'scalar','integer','positive'},'varName','objNumber'))
    addParameter(p,'errorCount',0, @(t) validateattributes(t, {'numeric'}, {'scalar','integer','positive'},'varName','errorCount'))
    addParameter(p,'verbose',1, @(t) validateattributes(t, {'numeric'}, {'scalar'},'varName','verbose'))

    parse(p,varargin{:})

    data_limits = p.Results.data_limits;
    r0set = p.Results.r0set;
    modifier = p.Results.modifier;
    useCentroid = logical(p.Results.useCentroid);
    objNumber = p.Results.objNumber;
    errorCount = p.Results.errorCount;
    verbose = logical(p.Results.verbose);

    % Using the centroid is only valid if the binned data is binary.
    % Silently turn useCentroid off if data is not binary.
    if ~islogical(binned_data) && any(binned_data(:)~=0 | binned_data(:)~=1)
        useCentroid = false;
    end

    if options.Use_Parallel && ~is_in_parallel()
        Use_Parallel = true;
    else
        Use_Parallel = false;
    end
end
