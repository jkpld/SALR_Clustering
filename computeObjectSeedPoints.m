function [seedPoints, Info] = computeObjectSeedPoints(binned_data, options, varargin)
% COMPUTEOBJECTSEEDPOINTS Compute the seed-points of an object.
%
% [seedPoints, Info] = computeObjectSeedPoints(binned_data, options)
% [seedPoints, Info] = computeObjectSeedPoints(binned_data, options, 'data_limits',dl,'r0set',r0,'modifier',m,'useCentroid',uc,'objNumber',on)
%
% Input parameters:
% binned_data : Binned data of a single object.
% options : An instance of class seedPointOptions.
%
% Optional parameters/value pairs:
% 'data_limits' : A 2xD array where the first row gives the minimum values 
%   of the data and the second row gives the maximum values of the data. D
%   is the dimension of the data. If 'data_limits' is not given, then the
%   size of each data bin will be assumed to be 1.
%       Ex. If the data is given by an array, dat, where each row is a
%       point and columns are dimensions, then the data_limits could be
%       computed using
%           data_limits = [min(dat,1); max(dat,1)];
% 'r0set' : NxD array of possible initial positions. Default is empty.
% 'modifier' : A matrix with size (size(binned_data) +
%   options.Potential_Padding_Size) that will be multiply the confining
%   potential before computing the confining force. If empty, then it is
%   ignored. Default is empty.
% 'useCentroid' : Logical flag. If true, then particles will not be
%   modeled, and the centroid of the binary mask will be returned as the
%   seed-point. (This is only useful in the context of locating nuclei
%   centers.) Default is false.
% 'objNumber' : The number of the current object. This is used with
%   computeNucleiCenters() for more usefull debug information.
%
% Output parameters:
% seedPoints : Nx2 array of computed seed points
% Info : A structure that always containes two fields, 
%   problem_scales : the problem_scales structure returned by
%     computeProblemScales(). 
%   message : integer giving an exit code
%     0 : everything is fine
%     1 : object is convex or too small (the center of the object
%       will be the seed point)
%     2 : less than 2 initial particles (the center of the object
%       will be the seed point)
%     3 : there was twice an error, object will be skipped
%
%   If options.debug is true, then Info will also have several other
%   fields:
%    r0 : initial points used in simulation
%    r_final : final location of simulated points
%    V : confining potential
%    ComputeInitialPointsInfo : structure of initial points information, 
%      see computeInitialPoints()


% See also MODELPARTICLEDYNAMICS EXTRACTCLUSTERCENTERS COMPUTEINITIALPOINTS
% CREATE_SCALEINVAR_CONFINING_POTENTIAL

% James Kapaldo

[binned_data, D, data_limits, r0set, modifier, useCentroid, objNumber, errorCount, Use_Parallel, verbose] = ...
    parse_inputs(binned_data, options, varargin{:});

Info = [];
DEBUG = options.Debug;

if useCentroid
    [seedPoints, Info] = computeCentroid(binned_data,data_limits,DEBUG,1);%'convex_or_tooSmall');
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
    Info.problem_scales = problem_scales; % save the problem scales.

    % Model dynamics -----------------------------------------------------
    % If there was only one initial particle, then we do not simulate it,
    % we will just return the centroid of the object.
    %
    % Note: Since we can model the object several times, if all of the
    % iterations have less than 4 particles and at least 1 of them only has
    % one particle, then we will not simulation it.
    num_r0 = cellfun(@(x) size(x,1), r0);
    if any(num_r0<2) && all(num_r0<4)
        [seedPoints,Info] = computeCentroid(binned_data,data_limits,DEBUG,2);%'lessThan_2_initial_particles');
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

    % Create a progress monitor
    progres = displayProgress(R, 10, verbose, Use_Parallel);
    Que = progres.start();

    if DEBUG
        r_final_n = cell(R,1);
        iteration_sizes = zeros(R,2);
        simulationInfo = cell(1,R);

        if Use_Parallel

            parfor n = 1:R
                [r_final_n{n}, simulationInfo{n}] = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r_final_n{n}, options);
                iteration_sizes(n,:) = [size(r_final_n{n},1), numel(cluster_sizes_n{n})];
                if ~isempty(Que), send(Que, obj), end
            end
        else
            for n = 1:R
                [r_final_n{n}, simulationInfo{n}] = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r_final_n{n}, options);
                iteration_sizes(n,:) = [size(r_final_n{n},1), numel(cluster_sizes_n{n})];
                progres.iteration_end()
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
                if ~isempty(Que), send(Que, obj), end
            end
        else
            for n = 1:R
                r = modelParticleDynamics(dV, r0{n}, options); % r_final is in solver space
                [seedPoints_n{n}, cluster_sizes_n{n}] = extractClusterCenters(r, options);
                progres.iteration_end()
            end
        end
    end

    delete(progres)

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
    end
    Info.message = 0;%'';

catch ME
    % It can be somtimes that there is an error in modelParticleDynamics
    % due to a particle flying out of the image region, or something that
    % rarely (hopefully) happens; therefore, we will try to run the
    % function again. If there is an error a second time, then we issue a
    % warning, save the error, and continue on to the next object.
    if any(strcmp('modelParticleDynamics',{ME.stack.name}))
        if errorCount < 1
            inputs.data_limits = data_limits;
            inputs.modifier = modifier;
            inputs.r0set = r0set;
            inputs.useCentroid = useCentroid;
            inputs.objNumber = objNumber;
            inputs.errorCount = errorCount + 1;
            [seedPoints, Info] = computeObjectSeedPoints(binned_data, options, inputs);
        else
            seedPoints = NaN(1,D);
            if DEBUG
                Info = emptyInfo(D);
            end
            Info.message = 3;
            Info.error = ME;
            fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', objNumber, objNumber)
            fprintf(2,'%s\n', getReport(Info.error,'basic'))
            fprintf('\n')
        end
    else
        rethrow(ME)
    end
end

end

function [seedPoint, Info] = computeCentroid(binned_data,data_limits,DEBUG,reason)

% Note, binned_data has not been padded and so the centroid is already in data
% units.
sz = size(binned_data);
D = numel(sz);
[BWpts{1:D}] = ind2sub(sz,find(binned_data));
BWpts = cat(2,BWpts{:});
seedPoint = mean(BWpts,1);

problem_scales = computeProblemScales([], sz, data_limits);
seedPoint = problem_scales.grid_to_data(seedPoint);

if DEBUG
    Info = emptyInfo(D);
end
Info.message = reason;
Info.problem_scales = problem_scales;

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
                  'message', NaN, ...
                  'problem_scales',struct());
end


function [binned_data, D, data_limits, r0set, modifier, useCentroid, objNumber, errorCount, Use_Parallel, verbose] = parse_inputs(binned_data, options, varargin)

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

    parse(p,varargin{:})

    data_limits = p.Results.data_limits;
    r0set = p.Results.r0set;
    modifier = p.Results.modifier;
    useCentroid = logical(p.Results.useCentroid);
    objNumber = p.Results.objNumber;
    errorCount = p.Results.errorCount;
    verbose = options.Verbose;

    % Using the centroid is only valid if the binned data is binary.
    % Silently turn useCentroid off if data is not binary.
    if ~islogical(binned_data) && any(binned_data(:)~=0 | binned_data(:)~=1)
        useCentroid = false;
    end

    currently_in_parallel = is_in_parallel();

    if options.Use_Parallel && ~currently_in_parallel
        Use_Parallel = true;
    else
        Use_Parallel = false;
    end

    % Get call stack to see if calling function is computNucleiCenters
    st = dbstack;
    if currently_in_parallel || any(contains({st.name},'nuclei','IgnoreCase',true))
        verbose = false;
    end
    
    if D == 2
        if ~isa(binned_data,'double')
            binned_data = double(binned_data);
        end
    end
end


%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
