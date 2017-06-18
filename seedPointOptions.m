classdef (ConstructOnLoad) seedPointOptions < matlab.mixin.CustomDisplay
    % SEEDPOINTOPTIONS  Set options needed for computing seed-point locations.
    %
    % Input can be structure array or parameter value pairs. Options not set
    % will be given default values. To see the default values, look at the
    % output with no inputs:
    %   defaultValues = seedPointOptions()
    %
    % seedPointOptions Properties:
    %
    % Point_Selection_Method - The method used to initialize the particle
    %   locations.
    %     'random' : Select N random locations from the binary mask where N
    %       is the area of the mask divided by the effective (hyper-)volume
    %       of a particle (pi*rs^2 in 2D), where rs is the
    %       Wigner_Seitz_Radius.
    %     'uniform' : A lattice is overlaid on the binary mask, and the
    %       centers of each lattice cell that are inside the binary mask
    %       are used as the initial points.
    %     'uniformRandom' : A lattice is overlaid on the binary mask, and
    %       from each lattice cell a random point from the mask is used as
    %       an initial point.
    %     'r0set_random' : From a set of possible initial positions, N
    %       random points are chosen, where N is the area of the mask
    %       divided by the effective volume of a particle (pi*rs^2 in 2D).
    %     'r0set_uniformRandom' : A lattice is overlaid on the binary mask,
    %       and from each lattice cell a random point from a set of
    %       possible initial positions is chosen.
    %   See computeInitialPoints() for more details.
    %   {'random', 'uniform', 'uniformRandom', 'r0set_random',
    %   'r0set_uniformRandom'}
    %
    % Wigner_Seitz_Radius - The effective size of each particle in the
    %   simulation. This sets the density of the particles used. The
    %   approximate number of particles used in a particular object will be
    %   object_size/particle_size.
    %   (0, Inf)
    %
    % Wigner_Seitz_Radius_Space - The coordinate space where the
    %   Wigner_Seitz_Radius is defined. For images of nuclei, 'data' space
    %   is the same thing as 'grid' space.
    %   {'data', 'grid', 'solver'}
    %
    % Maximum_Initial_Potential - The maximum allowed confining potential
    %   value for particles at their initial positions. Any possible
    %   initial position with a confining potential larger than this value
    %   will not be used as an initial position.
    %   (-Inf, Inf]
    %
    % Minimum_Initial_Potential - The minimum allowed confining potential
    %   value for particles at their initial positions. Any possible
    %   initial position with a confining potential smaller than this value
    %   will not be used as an initial position.
    %   [-Inf, Inf)
    %
    % Initial_Speed - The initial speed of the particles in the simulation.
    %   (-Inf, Inf)
    %
    %
    %
    % Potential_Parameters - [d0, r0, ra] The parameters describing the
    %   particle interaction potential: d0, the potential depth; r0, the
    %   location of the potential minimum; and ra, the attractive extent,
    %   which is the distance at which the potential goes from attractive
    %   to repulsive.
    %   [(-Inf,0], (0,Inf), (0,Inf) & > r0]
    %
    % Potential_Parameters_Space - The coordinate space where the
    %   Potential_Parameters are defined. (This value is not used.)
    %   {'data', 'solver'}
    %
    % Distance_Metric - The metric used for measuring distance between
    %   particles. If metric is Minkowski, then an extra agrument can be
    %   given with the exponent.
    %   {'euclidean'; 'cityblock'; 'chebychev'; {'minkowski', exponent}}
    %
    % Solver_Space_Attractive_Extent - The attractive extent of the
    %   particle interaction potential used when modeling the particle
    %   dynamics. If set to 'Attractive_Extent', then the attractive extent
    %   from the Potential_Parameters will be used. If set to a value, then
    %   the solver space will be isotropically scaled by a factor
    %   Solver_Space_Attractive_Extent/Potential_Parameters(3). Setting
    %   this parameter can give fine-tune control over the particle
    %   interaction.
    %   {'Attractive_Extent', (scalar, finite, real, positive)}
    %
    %
    %
    % Mass - The mass of the particles.
    %   (0, Inf)
    %
    % Coupling_Constant - The coupling constant value, k.
    %   (0, Inf)
    %
    % Charge_Normalization_Beta - The beta exponent for charge
    %   normalization based on number of particles.
    %   (-Inf, Inf)
    %
    %
    %
    % Potential_Type - The type of the confining potential used.
    %   {'distance_transform', 'density'}
    %
    % Potential_Modifier - A function handle that takes in the confining
    %   potential and outputs a modified potential with the same size. If
    %   empty, then it is ignored.
    %   {function_handle, []}
    %
    % Max_Distance_Transform - The maximum distance transform value. If
    %   set, then the object will be scaled so that the objects maximum
    %   distance transform is equal to this parameter. If NaN, or if
    %   Potential_Type is 'density', then this parameter is ignored.
    %   {(0, Inf), NaN}
    %
    % Max_Potential_Force - The maximum potential force. If set, then the
    %   potential force will be scaled so that its 99% value is equal to
    %   this parameter. If NaN, then it is ignored.
    %   {(0, Info), NaN}
    %
    % Potential_Padding_Size - The size of the padding to apply to an
    %   object mask. This is important to ensure the potential is defined
    %   some distance away from the object when running the simulation.
    %   (1, Inf)
    %
    %
    %
    % Iterations - The number of times to run the particle simulation. Each
    %   simulation will use a different set of initial positions. After
    %   running all iterations, the seed-points from each iteration will be
    %   clustered together. Using several iterations together with the
    %   Minimum_Cluster_Size can result in more stable/reproducable
    %   seed-points.
    %   [0, Inf), integer
    %
    % Minimum_Cluster_Size - The minimum number of particles in each
    %   particle cluster to be considered valid. The use of this parameter
    %   depends on the number of Iterations.
    %       If Iterations > 1, then the particles are the seed-points from
    %       each iteration. Therefore, any seed-point produced by less than
    %       Minimum_Cluster_Size iterations, will be removed. In this case,
    %       a good value of Minimum_Cluster_Size may be Iterations/3.
    %
    %       If Iterations = 1, then the particles are the particles of the
    %       the iteration. This means that each cluster of particles could
    %       potentially have N/C particles in them, if all particles were
    %       divided equally and where N is the number of particles in the
    %       simulation and C is the number of particle clusters
    %       (seed-points).
    %   [1, Inf), integer
    %
    % Particle_Damping_Rate - The rate at which the damping of the
    %   particles increases with simulation time.
    %   [0, Inf)
    %
    % Solver_Time_Range - Range over which the particles are simulated.
    %   This is passed to the ODE solver as the `tspan` parameter.
    %   [0, Inf)
    %
    % Maximum_Memory - The maximum amount of memory (per worker) that can
    %   be used to store the gradients of the confining potential. If the
    %   amount of space needed is above this limit, then a slower method
    %   will be used that does not take up as much memory.
    %   [0, Inf], [Gb]
    %
    %
    % Use_GPU - Determines if a GPU will be used to speed up calculation.
    %   logical
    %
    % Use_Parallel - Determines if multiple CPUs will be sued to speed up
    %   calculation.
    %   logical
    %
    % Verbose - Output information about progress of calculation.
    %   logical
    %
    % Debug - Output extra information about each iteration that may be
    %   helpful for debugging.
    %   logical


    % NOT USED ===============================================================
    % Curvature_Smoothing_Size - The standard deviation of the gaussians used
    %   for smoothing the boundary and calculating the curvature.
    %   [0, Inf) & integer , [pixels]
    %
    % Curvature_Max_Radius - Loosely the maximum object radius. The value will
    %   be used set a threshold for what is considered positive curvature
    %   (regions of positive curvature are concave regions of the boundary).
    %   The threshold will be set as 1/(2*Curvature_Max_Radius).
    %   [0, Inf) , [pixels]

    properties

        % Particle initialization
        
        Point_Selection_Method       = 'r0set_uniformRandom';
        Wigner_Seitz_Radius          = 5;
        Wigner_Seitz_Radius_Space    = 'grid'
        Maximum_Initial_Potential    = 1;
        Minimum_Initial_Potential    = -Inf;
        Initial_Speed                = 0.01;

        % Particle interaction
        
        Potential_Parameters         = [-1 2 15];
        Potential_Parameters_Space   = 'data'; 
        Distance_Metric              = 'euclidean';
        Solver_Space_Attractive_Extent = 'Attractive_Extent';

        % Particle parameters
        
        Mass                         = 1;
        Coupling_Constant            = 1;
        Charge_Normalization_Beta    = 1/3;

        % Confining potential parameters
        
        Potential_Type               = 'distance_transform'; 
        Potential_Modifier           = [];
        Max_Distance_Transform       = NaN;
        Max_Potential_Force          = NaN;
        Potential_Padding_Size       = 5;

        % Solver parameters
        
        Iterations                   = 1;
        Minimum_Cluster_Size         = 1;
        Particle_Damping_Rate        = 5e-4;
        Solver_Time_Range            = 0:10:1500;
        Maximum_Memory               = 1;

        % Computation options
        
        Use_Parallel                 = false;

        % Debug options
        
        Verbose                      = false;
        Debug                        = false;

    end

    properties (SetAccess = private)%, Hidden)
        ScaleInvarient_Potential_Extent = nan; % The problem will be scaled so that Potential_Extent is equal to this value and then it will be solved. This ensures the interaction potential has the same shape.
        ScaleInvarient_Potential_Minimum_Location

        potentialParameters % This holds a cache of all previously computed parameters
        potentialParameterIdx = [1 1 1];
        InteractionOptions = struct('type','SRALRR','params',[]);

        dist = 'euc';
        dist_arg = [];
    end

    properties (Hidden)
        Use_GPU        = false;
        Use_ConvexHull = true;
    end

    methods
        function options = seedPointOptions(varargin)
            
            % First get the potential paramters
            try
                options.potentialParameters = load('potentialParameters.mat');
            catch
                error('seedPointOptions:missingPotentialParamters','The file potentialParameters.mat is missing or not on the path. This file is required for setting the potentialParamters. It may be created with the function computePotentialParamters if you do not have it.')
            end

            % Make sure the potential parameters are initialized correctly
            options.Solver_Space_Attractive_Extent = options.Solver_Space_Attractive_Extent;
            options = setPotentialParameter(options, options.Potential_Parameters);

            % Now assign any properties given.
            if nargin > 0
                if numel(varargin) == 1 && isstruct(varargin{1})
                    op = varargin{1};
                    fd = fieldnames(op)';
                    for fieldname = fd
                        options.(char(fieldname)) = op.(char(fieldname));
                    end
                else
                    if ~mod(nargin+1,2)
                        error('seedPointOptions:badInput','Input must be structure with properties as fields, or property/value list with even number of elements.')
                    end
                    for i = 1:2:numel(varargin)
                        options.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end

        function obj = set.Wigner_Seitz_Radius(obj,value)
            validateattributes(value,{'double'},{'positive','scalar','real','finite'})
            obj.Wigner_Seitz_Radius = value;
        end

        function obj = set.Wigner_Seitz_Radius_Space(obj,value)

            methods = {'data'; 'grid'; 'solver'};

            i = find(strncmpi(value,methods,length(value)));
            if length(i) > 1
                error('seedPointOptions:AmbiguousWignerSeitzRadiusSpace', 'Ambiguous space, %s, valid options are ''data'', ''grid'', and ''solver''.', value);
            elseif isempty(i)
                error('seedPointOptions:UnknownWignerSeitzRadiusSpace', 'Invalid space, %s. Valid options are ''data'', ''grid'', and ''solver''.', value);
            end

            obj.Wigner_Seitz_Radius_Space = lower(methods{i});
        end

        function obj = set.Initial_Speed(obj,value)
            validateattributes(value,{'double'},{'scalar','real','finite'})
            obj.Initial_Speed = value;
        end

        function obj = set.Mass(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            obj.Mass = value;
        end

        function obj = set.Coupling_Constant(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            obj.Coupling_Constant = value;
        end

        function obj = set.Point_Selection_Method(obj,value)
            value = validatestring(value,{'random','uniform','uniformRandom','r0set_random','r0set_uniformRandom'});
            obj.Point_Selection_Method = value;
        end

        function obj = set.Potential_Parameters(obj, value)
            validateattributes(value,{'double'},{'numel',3,'real','finite'})

            if value(1) > 0
                error('seedPointOptions:badSet','The potential depth must less than zero.')
            end
            if value(3) <= value(2)
                error('seedPointOptions:badSet','The potential minimum location must be smaller than the potential extent.')
            end

            obj.Potential_Parameters = value;

            if ischar(obj.Solver_Space_Attractive_Extent) && ...
                    strcmp(obj.Solver_Space_Attractive_Extent,'Attractive_Extent') %#ok<MCSUP>
                obj.ScaleInvarient_Potential_Extent = value(3); %#ok<MCSUP>
            end

            obj = setPotentialParameter(obj, value);
        end

        function obj = set.Potential_Parameters_Space(obj,value)

            validOptions = {'data', 'solver'};
            idx = find(strncmpi(value, validOptions, length(value)));

            if isempty(idx)
                error('seedPointOptions:unknownPtntlPrmtrSpace','The Potential_Parameters_Space, %s, is not valid. Valid options are ''data'' and ''solver''.',value)
            elseif length(idx)>1
                error('seedPointOptions:unknownPtntlPrmtrSpace','The Potential_Parameters_Space, %s, is ambiguous. Please be more specific.',value)
            end

            obj.Potential_Parameters_Space = validOptions{idx};

            if idx == 2
                obj.Solver_Space_Attractive_Extent = 'Attractive_Extent'; %#ok<MCSUP>
            end
        end

        function obj = set.Potential_Type(obj, value)

            validOptions = {'distance_transform','density'};
            idx = find(strncmpi(value, validOptions, length(value)));

            if isempty(idx)
                error('seedPointOptions:unknownPtntlType','The Potential_Type, %s, is not valid. Valid options are ''distance_transform'' and ''density''.',value)
            elseif length(idx)>1
                error('seedPointOptions:unknownPtntlType','The Potential_Type, %s, is ambiguous. Please be more specific.',value)
            end

            obj.Potential_Type = validOptions{idx};
        end

        function obj = set.Potential_Modifier(obj, value)
            if isempty(value) || isa(value,'function_handle')
                obj.Potential_Modifier = value;
            else
                error('seedPointOptions:invalidModifier','The potential modifier must either be a function handle or empty.')
            end
        end

        function obj = set.Max_Distance_Transform(obj, value)

            validateattributes(value,{'double'},{'scalar','positive','real'})
            if isnan(value)
                if strcmp(obj.Potential_Parameters_Space,'max_distance_transform') %#ok<MCSUP>
                    warning('seedPointOptions:nanMaxDT','Setting Potential_Parameters_Space to ''data''. Potential_Parameters_Space cannot be set to ''max_distance_transform'' when Max_Distance_Transform is NaN.')
                    obj.Potential_Parameters_Space = 'data'; %#ok<MCSUP>
                end
            elseif isinf(value)
                error('seedPointOptions:infMaxDT','Max_Distance_Transform must be a scalar positive real finite number, or NaN.')
            end

            obj.Max_Distance_Transform = value;
        end

        function obj = set.Max_Potential_Force(obj, value)

            validateattributes(value,{'double'},{'scalar','positive','real'})
            if isinf(value)
                error('seedPointOptions:infMaxForce','Max_Potential_Force must be a scalar positive real finite number, or NaN.')
            end

            obj.Max_Potential_Force = value;
        end

        function obj = set.Potential_Padding_Size(obj,value)
            validateattributes(value,{'double'},{'integer','scalar','nonnegative','real','finite'})
            obj.Potential_Padding_Size = value;
        end

        function obj = set.Iterations(obj, value)
            validateattributes(value, {'numeric'}, {'integer','positive','real','finite'})
            obj.Iterations = value;
        end

        function obj = set.Minimum_Cluster_Size(obj, value)
            validateattributes(value, {'double'}, {'scalar','nonnegative','finite','real'})
            obj.Minimum_Cluster_Size = value;
        end

        function obj = set.Maximum_Initial_Potential(obj,value)
            validateattributes(value,{'double'},{'scalar','real','>',-Inf,'<=',Inf})
            if value < obj.Minimum_Initial_Potential %#ok<MCSUP>
                error('seedPointOptions:badInput','MaximumM_Initial_Potential must be larger than inimum_Initial_Potential.')
            end
            obj.Maximum_Initial_Potential = value;
        end

        function obj = set.Minimum_Initial_Potential(obj,value)
            validateattributes(value,{'double'},{'scalar','real','>=',-Inf,'<',Inf})
            if value > obj.Maximum_Initial_Potential %#ok<MCSUP>
                error('seedPointOptions:badInput','Minimum_Initial_Potential must be smaller than Maximum_Initial_Potential.')
            end
            obj.Minimum_Initial_Potential = value;
        end

        function obj = set.Distance_Metric(obj, input)

            if iscell(input)
                if length(input) > 1
                    extraArg = input{2};
                else
                    extraArg = [];
                end
                metric = input{1};
            else
                metric = input;
                extraArg = [];
            end

            additionalArg = [];

            % The code below mostly comes from Matlab's pdist function.
            methods = {'euclidean'; 'cityblock'; 'chebychev'; 'minkowski'};

            i = find(strncmpi(metric,methods,length(metric)));
            if length(i) > 1
                error('seedPointOptions:AmbiguousDistance', 'Ambiguous distance, %s, try entering the full distance metric''s name.', metric);
            elseif isempty(i)
                error('seedPointOptions:UnknownDistance', 'Unknown distance metric, %s.', metric);
            else
                metric = lower(methods{i}(1:3));
            end

            if strcmp(metric,'min') % Minkowski distance
                additionalArg = extraArg;
                if isempty(additionalArg)
                    metric = 'euc';
                    i = 1;
                    additionalArg = [];
                elseif ~( isscalar(additionalArg) && additionalArg > 0)
                    error('seedPointOptions:InvalidExponent','Invalid exponent for Minkowski metric.');
                elseif isinf(additionalArg) %the exponent is inf
                    metric = 'che';
                    i = 3;
                    additionalArg = [];
                elseif additionalArg == 2 %the exponent is 2
                    metric = 'euc';
                    i = 1;
                    additionalArg = [];
                elseif additionalArg == 1 %the exponent is 1
                    metric = 'cit';
                    i = 2;
                    additionalArg = [];
                end
            end

            obj.dist = metric; %#ok<MCSUP>
            obj.dist_arg = additionalArg; %#ok<MCSUP>

            if isempty(additionalArg)
                obj.Distance_Metric = methods{i};
            else
                obj.Distance_Metric = {methods{i}, additionalArg};
            end

        end

        function obj = set.Use_GPU(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('seedPointOptions:badInput','Expected input to be logical.')
            end
            obj.Use_GPU = value;
        end

        function obj = set.Use_Parallel(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('seedPointOptions:badInput','Expected input to be logical.')
            end

            if value == 1
                % Make sure the parallel computing toolbox is installed and
                % can be used.
                ver_info = ver('distcomp');
                ver_lic = license('test','Distrib_Computing_Toolbox');
                if isempty(ver_info) || ~ver_lic
                    warning('seedPointOptions:NoParallelToolbox','Cannot use parallel computing. The parallel computing toolbox either is not installed or there is no license for it.')
                    obj.Use_Parallel = false;
                else
                    obj.Use_Parallel = true;
                end
            else
                obj.Use_Parallel = false;
            end
        end

        function obj = set.Maximum_Memory(obj,value)
            validateattributes(value,{'double'},{'positive','scalar','real'})
            obj.Maximum_Memory = value;
        end

        function obj = set.Verbose(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('seedPointOptions:badInput','Expected input to be logical.')
            end
            obj.Verbose = value;
        end

        function obj = set.Debug(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('seedPointOptions:badInput','Expected input to be logical.')
            end
            obj.Debug = value;
        end

        function obj = set.Solver_Time_Range(obj,value)
            validateattributes(value,{'double'},{'nonnegative','real','finite'})
            obj.Solver_Time_Range = value;
        end

        function obj = set.Particle_Damping_Rate(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','real','finite'})
            obj.Particle_Damping_Rate = value;
        end

        function obj = set.Charge_Normalization_Beta(obj,value)
            validateattributes(value,{'double'},{'scalar','real','finite'})
            obj.Charge_Normalization_Beta = value;
        end

        function obj = set.Solver_Space_Attractive_Extent(obj, value)
            if ischar(value)
                if ~strcmp(value, 'Attractive_Extent')
                    error('seedPointOptions:UnknownInput', 'Unknown input, %s, for Solver_Space_Attractive_Extent. Value should be the string "Attractive_Extent", or a scalar, real, positive, finite value.', value);
                else
                    obj.Solver_Space_Attractive_Extent = value;
                    obj.ScaleInvarient_Potential_Extent = obj.Potential_Parameters(3); %#ok<MCSUP>
                    obj = setPotentialParameter(obj,obj.Potential_Parameters); %#ok<MCSUP>
                end
            else
                validateattributes(value,{'double'},{'scalar','real','finite','positive'})
                if strcmp(obj.Potential_Parameters_Space, 'solver') %#ok<MCSUP>
                    warning('seedPointOptions:overspecifiedParameters', 'Solver space will be overspecified by setting Solver_Space_Attractive_Extent since Potential_Parameters_Space is set to ''solver''. If you want to set Solver_Space_Attractive_Extent to a value other than the Potential parameters attractive extent, then first change the Potential_Parameters_Space to ''data''.')
                else
                    obj.Solver_Space_Attractive_Extent = value;
                    obj.ScaleInvarient_Potential_Extent = value; %#ok<MCSUP>
                    obj = setPotentialParameter(obj,obj.Potential_Parameters); %#ok<MCSUP>
                end
            end
        end

        function out = plotPotential(obj,maxR)

            if nargin < 2
                maxR = 1.3;
            end

            % Compute the interaction potential
            scaleFactor = obj.Potential_Parameters(3) / obj.ScaleInvarient_Potential_Extent;

            x = obj.InteractionOptions.params;
            r = 0:0.05:maxR*obj.ScaleInvarient_Potential_Extent;
            Vint = 1./(r+0.2) - x(1)*exp(-(r-x(2)).^2/(2*x(3)^2));
            r = r * scaleFactor;

            if nargout > 0
                out = {r,Vint, r(1:end-1)+0.025,diff(Vint)};
                return
            end

            figure

            % Plot interaction potential
            subplot(2,1,1)
            line([r(1),r(end)], [0 0],'linestyle','--','color','k')
            line(r,Vint,'color','b','linewidth',2)

            xticks = [0, obj.Potential_Parameters(2), obj.Potential_Parameters(3)];
            set(gca,'XTick',xticks,'XLim',[r(1),r(end)],'YLim',1.2*abs(obj.Potential_Parameters(1))*[-1,1]);
            title(sprintf('Interaction potential\n (d_0=%0.2f, r_0=%0.2f, r_a=%0.2f) @ r_{a,SI}=%0.2f', obj.Potential_Parameters(1), obj.Potential_Parameters(2),obj.Potential_Parameters(3), obj.ScaleInvarient_Potential_Extent))

            % Plot interaction force
            subplot(2,1,2)
            line([r(1),r(end)], [0 0],'linestyle','--','color','k')
            line(r(1:end-1)+0.025,diff(Vint),'color','b','linewidth',2)

            set(gca,'XTick',xticks,'XLim',[r(1),r(end)],'YLim',1.2*max(diff(Vint))*[-1,1]);
            title('Interaction force')

            % Set theme
            setTheme(gcf,'light')
        end

        function validateInteractionPotential(obj)
            % Confirm that the attractive extent, depth, and center
            % actually are where we set them to be. Not all combinations of
            % parameters are possible (small center and large extent for
            % example)

            x = obj.InteractionOptions.params;
            Vint = @(r) 1./(r+0.2) - x(1)*exp(-(r-x(2)).^2/(2*x(3)^2));
            dVint = @(D) -1./(D + 0.2).^2 + (x(1)*(D-x(2))/(x(3)^2)) .* exp(-(D-x(2)).^2/(2*x(3)^2));

            % Test extent
            y = obj.ScaleInvarient_Potential_Extent;
            y_hat = fzero(dVint, y*2);
            y_err = abs(y_hat - y)/y;

            if y_err > 1e-2
                warning('seedPointOptions:unsolvableInteractionPotential', 'The potential attractive extent is %0.2f%% different than the set potential attractive extent. Consider modifying the potential parameters to ensure the interaction potential is as expected.',y_err*100)
            end

            % Test minimum location
            y = obj.ScaleInvarient_Potential_Minimum_Location;
            y_hat = fzero(dVint, y);
            y_err = abs(y_hat - y)/y;

            if y_err > 1e-2
                warning('seedPointOptions:unsolvableInteractionPotential', 'The potential minimum location is %0.2f%% different than the set potential minimum location. Consider modifying the potential parameters to ensure the interaction potential is as expected.',y_err*100)
            end

            % Test depth
            y = obj.Potential_Parameters(1);
            y_hat = Vint(y_hat);
            y_err = abs(y_hat - y)/y;

            if y_err > 1e-2
                warning('seedPointOptions:unsolvableInteractionPotential', 'The potential depth is %0.2f%% different than the set potential depth. Consider modifying the potential parameters to ensure the interaction potential is as expected.',y_err*100)
            end
        end
    end

    methods (Access = private)
        function obj = setPotentialParameter(obj,params)

            depth = params(1);
            center = params(2);
            extent = params(3);

            center = obj.ScaleInvarient_Potential_Extent * center / extent;
            obj.ScaleInvarient_Potential_Minimum_Location = center;

            extent = obj.ScaleInvarient_Potential_Extent;

            depth_idx = find(obj.potentialParameters.depth == depth);
            center_idx = find(obj.potentialParameters.center == center);
            extent_idx = find(obj.potentialParameters.extent == extent);

            if isempty(depth_idx) || isempty(center_idx) || isempty(extent_idx)
                % Need to compute a new set of parameters

                parameterStruct = computePotentialParameters(depth,center,extent);

                if isempty(depth_idx)
                    obj.potentialParameters.depth = [obj.potentialParameters.depth, depth];
                    depth_idx = numel(obj.potentialParameters.depth);
                end
                if isempty(center_idx)
                    obj.potentialParameters.center = [obj.potentialParameters.center, center];
                    center_idx = numel(obj.potentialParameters.center);
                end
                if isempty(extent_idx)
                    obj.potentialParameters.extent = [obj.potentialParameters.extent, extent];
                    extent_idx = numel(obj.potentialParameters.extent);
                end

                obj.potentialParameters.parameters = [obj.potentialParameters.parameters; parameterStruct.parameters];
                obj.potentialParameters.parametersIdx = [obj.potentialParameters.parametersIdx; encodePotentialIdx(depth_idx, center_idx, extent_idx)];
                obj.InteractionOptions.params = parameterStruct.parameters;
            else

                idx = obj.potentialParameters.parametersIdx == encodePotentialIdx(depth_idx, center_idx, extent_idx);

                if ~any(idx)
                    % Do not have this combination of parameters, need to
                    % compute them
                    parameterStruct = computePotentialParameters(depth,center,extent);

                    obj.potentialParameters.parameters = [obj.potentialParameters.parameters; parameterStruct.parameters];
                    obj.potentialParameters.parametersIdx = [obj.potentialParameters.parametersIdx; encodePotentialIdx(depth_idx, center_idx, extent_idx)];
                    obj.InteractionOptions.params = parameterStruct.parameters;
                else
                    % Already have the parameters, just need to set them
                    obj.InteractionOptions.params = obj.potentialParameters.parameters(idx,:);
                end
            end

            obj.potentialParameterIdx = [depth_idx, center_idx, extent_idx];

            validateInteractionPotential(obj)
        end
    end

    methods (Access = protected)
        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            fprintf('\n%s\n',[className,' with properties:']);

            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
            %          for i = 1:length(propgroup)
            %              propgroup(i).Aligned = 0;
            %              matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup(i))
            %              fprintf('\n')
            %          end


        end

        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                % property groups for scalars
                nGroups = 8;
                %             gTitles = {'Particle initialization';
                %                        'Particle interaction';
                %                        'Particle parameters';
                %                        'Confining potential parameters';
                %                        'Solver parameters';
                %                        '2D mask clean up';
                %                        'Computation options';
                %                        'Debug options'};
                gTitles = repmat({' '},1,nGroups);

                propLists = cell(1,nGroups);
                propLists{1} = {'Point_Selection_Method';
                    'Wigner_Seitz_Radius';
                    'Maximum_Initial_Potential';
                    'Minimum_Initial_Potential';
                    'Initial_Speed'};
                propLists{2} = {'Potential_Parameters';
                    'Potential_Parameters_Space';
                    'Distance_Metric';
                    'Solver_Space_Attractive_Extent'};
                propLists{3} = {'Mass';
                    'Coupling_Constant';
                    'Distance_Metric'};
                propLists{4} = {'Max_Distance_Transform';
                    'Max_Potential_Force';
                    'Potential_Padding_Size'};
                propLists{5} = {'Particle_Damping_Rate';
                    'Solver_Time_Range';
                    'Maximum_Memory'};
                propLists{6} = {'Minimum_Hole_Size'};
                propLists{7} = {'Use_GPU';
                    'Use_Parallel'};
                propLists{8} = {'Debug';
                    'Object_Of_Interest'};


                for i = nGroups:-1:1
                    propgrp(i) = matlab.mixin.util.PropertyGroup(propLists{i},gTitles{i});
                end
            end
        end
    end
end

function out = encodePotentialIdx(depth_idx,center_idx,extent_idx)
out = uint32(depth_idx + bitshift(center_idx,8) + bitshift(extent_idx,16));
end
