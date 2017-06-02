function [seedPoints,Info] = compute_seedPoints_nd(BW,data_range,options,varargin)



[data_type, r0_set, data_range, N, minClusterSize, verbose] = parse_inputs(BW, data_range, varargin{:});

PAD_SIZE = options.Potential_Padding_Size;

D = ndims(BW);
sz = size(BW);
grid_range = sz - 1;



Use_Parallel = options.Use_Parallel;
if N < 2
    Use_Parallel = false;
end


% Create confining potential ---------------------------------------------
switch data_type
    case 'density'
        V = min(1.5,1./BW);
        V = padarray(V, PAD_SIZE*ones(1,D), 1.5);
%         BW = V < 0.9;
%         dt = bwdist(BW);
%         V(~BW) = dt(~BW);
%         clear dt
    case 'mask'
        [V, scale_factor] = create_scaleInvar_confining_potential(BW,options);
        options.Scale_Factor = scale_factor;
end

% Compute scale factors --------------------------------------------------
[grid_spacing, grid_to_solver, grid_to_data] = computeGridScaleFactors(options, grid_range, data_range);
solver_to_data = @(x) (x./grid_to_solver - PAD_SIZE) .* grid_to_data;
data_to_grid = @(x) x./grid_to_data;

% Compute confining force ------------------------------------------------
dV = confiningForce(V, grid_spacing, grid_to_solver, options);

% Get initial positions --------------------------------------------------

BW = (V < options.Maximum_Initial_Potential) & (V > options.Minimum_Initial_Potential);
rsGrid = getGridWignerSeitz(grid_to_data, grid_to_solver, options);

if isempty(r0_set)
    r0_set = computeInitialPoints(BW, options.Point_Selection_Method, rsGrid, [] ,N);
else
    % Pad the initial posiitons if the were given
    r0_set = r0_set + options.Potential_Padding_Size;
    r0_set = repmat({r0_set},1,N);
end

% Convert initial positions to solver space
r0_set = cellfun(@(x) x .* grid_to_solver, r0_set, 'UniformOutput',false);

% Model dynamics ---------------------------------------------------------
seedPoints_n = cell(N,1);
Info_n = cell(N,1);


if Use_Parallel
    warning('off','MATLAB:structOnObject')
    options = struct(options);
    warning('on','MATLAB:structOnObject')

    parfor n = 1:N
        [r, Info_n{n}] = modelParticleDynamics(dV,r0_set{n},options); % r is in solver space
        [seedPoints_n{n}, clstSz] = extractClusterCenters(r,options);

        Info_n{n}.r0 = solver_to_data(r0_set{n});
        Info_n{n}.r_final = solver_to_data(r);
        Info_n{n}.seedPoints = solver_to_data(seedPoints_n{n});
        Info_n{n}.cluster_size = clstSz;
    end
else
    if verbose
        generateDisplayAt = unique(round(linspace(1,N,7)));
    else
        generateDisplayAt = nan;
    end
    processTimes = zeros(1,N);

    for n = 1:N
        procTime = tic;

        [r, Info_n{n}] = modelParticleDynamics(dV,r0_set{n},options); % r is in solver space
        [seedPoints_n{n}, clstSz] = extractClusterCenters(r,options);

        Info_n{n}.r0 = solver_to_data(r0_set{n});
        Info_n{n}.r_final = solver_to_data(r);
        Info_n{n}.seedPoints = solver_to_data(seedPoints_n{n});
        Info_n{n}.cluster_size = clstSz;

        processTimes(n) = toc(procTime);
        if any(n == generateDisplayAt)
            fprintf('%s >> Iteration %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),n,N,sum(processTimes(1:n))/60,mean(processTimes(1:n))*N/60)
        end
    end
end

% Combine iterations and remove clusters with less than minClusterSize
% particles. -------------------------------------------------------------
seedPoint_set = cat(1,seedPoints_n{:}); % This is in solver space

% Remove any particles at a confining potential larger than 1, these
% particles are outside of the object - this can occure when using density
% as the confining potential because we do not have a quadratic confining
% potential around the object like we do with the distance transform
% confining potential.
Vi = interpolatePotential(V,seedPoint_set,grid_to_solver);
seedPoint_set(Vi>1,:) = [];


if isempty(seedPoint_set)
    error('compute_seedPoints_nd:strangeError','All computed seed-points lie outside of the object (V>1).')
end

if N > 1
    [seedPoints, clstSz] = extractClusterCenters(seedPoint_set, options);
else
    seedPoints = seedPoint_set;
    clstSz = Info_n{1}.cluster_size;
end
toRemove = clstSz < minClusterSize;
seedPoints(toRemove,:) = [];

seedPoints = solver_to_data(seedPoints);

% Create Info structure for output ---------------------------------------
Info = [];
Info.iteration_info = Info_n;
Info.seedPoint_set = solver_to_data(seedPoint_set);
Info.cluster_sizes = clstSz;
Info.solver_to_data = solver_to_data;
Info.data_to_grid = data_to_grid;
end

function [dataType, r0_set, data_range, iterations, minClusterSize, verbose] = parse_inputs(inputData,data_range,varargin)

% Get input values
p = inputParser;
p.FunctionName = 'compute_seedPoints_nd';

if any(inputData(:)<0)
    error('compute_seedPoints_nd:negative','The input data must be non-negative.')
end

validateattributes(data_range,{'numeric'},{'numel',ndims(inputData),'nonnegative','finite','real'})

inputClass = class(inputData);
switch inputClass
    case 'logical'
        defaultType = 'mask';
    case {'single','double'}
        defaultType = 'density';
    otherwise
        error('compute_seedPoints_nd:unknownInputType','The input data should be logical or float (single, double).')
end

addParameter(p,'iterations',1, @(t) validateattributes(t,{'numeric'},{'integer','scalar','positive','real','finite'}));
addParameter(p,'dataType',defaultType, @(t) any(strcmpi(t,{'','mask','density'})))
addParameter(p,'minClusterSize',0, @(t) validateattributes(t,{'numeric'},{'scalar','nonnegative','real','finite'}));
addParameter(p,'r0_set',[], @(t) isempty(t) || validateattributes(t,{'numeric'},{'2d','real','finite'}))
addParameter(p,'verbose',1, @(t) validateattributes(t,{'numeric'},{'scalar'}))

parse(p,varargin{:})

iterations = p.Results.iterations;
dataType = p.Results.dataType;
minClusterSize = p.Results.minClusterSize;
r0_set = p.Results.r0_set;
verbose = logical(p.Results.verbose);

end

% function r0_set = computeInitialPoints_nd(BW,options,N)
% 
% if nargin < 3
%     N = 1;
% end
% sz = size(BW);
% D = length(sz);
% 
% % The lattice constant is the volume of a D-dimension ball to the (1/D)
% % power:
% a = pi^(1/2)/gamma(D/2+1)^(1/D) * getGridWignerSeitz(grid_to_data, grid_to_solver, options);
% 
% switch options.Point_Selection_Method
% 
%     case 'random'
%         idx = find(BW);
% 
%         M = round(numel(idx)/a^D); % Number of particles to use
% 
%         r0_set = cell(1,N);
%         for i = 1:N
%             r0_set_i = cell(1,D);
%             idx = idx(randperm(size(idx,1),M));
%             [r0_set_i{:}] = ind2sub(sz,idx);
%             r0_set{i} = cat(2,r0_set_i{:});
%         end
% 
%     case 'uniform'
%         r0_set = cell(1,D);
% 
%         for i = 1:D
%             r0_set{i} = 1:a:sz(i);
%         end
% 
%         [r0_set{:}] = ndgrid(r0_set{:});
%         r0_set = cellfun(@(x) x(:), r0_set,'UniformOutput',false);
% 
%         ind = sub2ind(sz,round(r0_set{:}));
%         remove = BW(ind) == 0
% 
%         r0_set = cat(2,r0_set{:});
%         r0_set(remove,:) = [];
%         r0_set = repmat({r0_set},1,N);
% 
%     case 'uniformRandom'
% 
%         idx = find(BW);
%         locs = cell(1,D);
%         [locs{:}] = ind2sub(sz,idx);
% 
%         bin = zeros(numel(idx),D,'int8');
%         for i = D:-1:1
%             edges = 0 : a : (sz(i)+a);
%             bin(:,i) = discretize(locs{i},edges);
%         end
% 
%         locs = cat(2,locs{:});
% 
%         r0_set = cell(1,N);
%         for i = 1:N
%             r0_set_idx = accumarray(bin(all(bin>0,2),:), (1:numel(idx))', [], @(x) x(randi(numel(x),1)), 0);
%             r0_set_idx = r0_set_idx(logical(r0_set_idx));
%             r0_set{i} = locs(r0_set_idx,:);
%         end
%     otherwise
%         error('compute_seedPoints_nd:notValidSelectionMethod','The point selection method, %s, is either unknown or has not been programed from N-D data. Valid options are ''random'', ''uniform'', and ''uniformRandom''.',options.Point_Selection_Method)
% end
% 
% end

function Vi = interpolatePotential(V,r,grid_to_solver)

sz = size(V);
D = size(r,2);
dx = cell(1,D);

for i = 1:D
    dx{i} = (1:sz(i)) * grid_to_solver(i);
end

Vi = griddedInterpolant(dx,V);
Vi = Vi(r);
end
