classdef seedPointOptions
% SEEDPOINTOPTIONS  Set options needed for computing seed point locations.
%
% Input can be structure array or parameter value pairs. Options not set
% will be given default values. To see the default values, look at the
% output with no inputs:
%   defaultValues = seedPointOptions()
%
% seedPointOptions Properties:
%
% Solver_Space_Scale - The scale defining the space the particles
%   are modeled in. This is the space where the interaction potential
%   parameters (Extent, Minimum Location, and Depth) are defined. There are
%   four options:
%       'grid' : The problem is directly solved in the data-storage space
%       'data' : The problem is solved in the data space
%       'distance_transform' : The data is solved in 'grid' space, if the
%           distance transform is not scaled (Max_Distance_Transform =
%           NaN). If the distance transform is scaled
%           (Max_Distance_Transform = 'scalar'), then the 'grid' space is
%           scaled by the same ratio as the distance transform. - This is
%           the spaced used for the scale invariant seed-point calculations
%           of the nuclei centers.
%       [~>10, Inf) : Directly set a scalar scale factor
%   {[~>10, Inf), 'distance_transform', 'grid', 'data'}
%           # Comment : The interaction potential cannot be solved well
%           when the Extent and Minimum Location are small (Extent ~< 5 and
%           Minimum Location ~< 1); this is why the Solver_Space_Scale
%           should be >~10
%
% Data_Storage_Centers - A cell array of arrays giving the centers of each
%   data bin in data units.
%
% Wigner_Seitz_Radius - The effective size of each particle in the
%   simulation. This sets the density of the particles used. The
%   approximate number of particles used in a particular object will be
%   area(object)/(pi* r_s^2), where r_s is the Wigner_Seitz_Radius.
%   (0, Inf)
%
% Potential_Depth - The depth of the potential.
%   (-Inf,0]
%
% Potential_Minimum_Location - The location of the potential minimum.
%   (0, Inf)
%
% Potential_Extent - The radius at which the potential goes from attractive
%   to repulsive.
%   (0, Inf) & > Potential_Minimum_Location
%
% Potential_Padding_Size - The size of the padding to apply to an object
%   mask. This is important so that the potential is defined some distance
%   away from the object when running the simulation.
%   (1, Inf)
%
% Maximum_Initial_Potential - The maximum confining potential value that a
%   particle can have as its initial position. Any possible initial
%   position with a confining potential larger than this value will not be
%   used as an initial particle position.
%   (-Inf, Inf]
%
% Minimum_Initial_Potential - The minimum confining potential value that a
%   particle can have as its initial position. Any possible initial
%   position with a confining potential smaller than this value will not be
%   used as an initial particle position.
%   (-Inf, Inf]
%
% Potential_Scale - Value used to set the depth of the confining potential
%   well. Set to NaN to use the natural depth.
%   {(0, Inf), NaN}
%
% Distance_Metric - The metric used for measuring distance between
%   particles. If metric is Minkowski, then an extra agrument can be given
%   with the exponent.
%   {'euclidean'; 'cityblock'; 'chebychev'; {'minkowski', exponent}}
%
% Initial_Speed - The initial speed of the particles in the simulation.
%   (-Inf, Inf)
%
% Mass - The mass of the particles.
%   (0, Inf)
%
% Coupling_Constant - The coupling constant value, k.
%   (0, Inf) 
%
% Particle_Damping_Rate - The rate at which the damping of the particles
%   increases with simulation time.
%   [0, Inf)
%
% Charge_Normalization_Beta - The beta exponent for charge normalization
%   based on number of particles.
%   (-Inf, Inf)
%
% Solver_Time_Range - Range over which the particles are simulated.
%   [0, Inf) 
%
% Point_Selection_Method - The method used to initialize the particle
%   locations.
%   {'random','uniform','uniformRandom','r0set_random','r0set_uniformRandom'}
%
% Minimum_Hole_Size - The minimum hole size (area) allowed in the mask.
%   [0, Inf)
%
% Use_GPU - Determines if a GPU will be used to speed up calculation.
%   logical
%
% Use_Parallel - Determines if multiple CPUs will be sued to speed up
%   calculation.
%   logical
%
% Maximum_Memory - The maximum amount of memory that can be used to store
%   the gradients of the confining potential. The of the amount of space
%   needed is above this limit, then a slower method will be used that does
%   not take up memory.
%   [0, Inf], [Gb]
%
% Debug -
%   Determines if additional information will be returned from each
%   function.
%   logical
%
% Object_Of_Interest - The index of an object of interest to declump. If an
%   object is giving an error. Then set this property to the object index
%   and set Debug to true. Plots showing intermediate steps of the
%   calculation will be shown that should help you debug the problem. Leave
%   this property empty for normal use.


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

        Solver_Space_Scale           = 'auto';
        Data_Storage_Centers
        
        % Initial particles
        Wigner_Seitz_Radius          = 10;
        Point_Selection_Method       = 'r0set_uniformRandom';
        Maximum_Initial_Potential    = Inf;
        Minimum_Initial_Potential    = -Inf;
        
        % Particle interaction
        Potential_Depth
        Potential_Minimum_Location
        Potential_Extent
        Potential_Parameter_Space    = 'data'; % {'data', 'max_distance_transform'}
        Distance_Metric              = 'euclidean';
        
        % Particle parameters
        Initial_Speed                = 0.01;
        Mass                         = 1;
        Coupling_Constant            = 1;
        Particle_Damping_Rate        = 5e-4;
        Charge_Normalization_Beta    = 1/3;
        
        % Potential parameters
        Potential_Padding_Size       = 5;
        Potential_Scale              = NaN;
        
        % Solver parameters
        Solver_Time_Range            = 0:10:1500;
        Maximum_Memory               = 1;
        
        % 2D mask clean up
        Minimum_Hole_Size            = 50;

        % Computation options
        Use_GPU                      = false;
        Use_Parallel                 = false;
        
        % Debug options
        Debug                        = false;
        Object_Of_Interest           = [];
        
        % Scale invariant interaction parameters
        ScaleInvarient_Potential_Extent = 15; % The problem will be scaled so that Potential_Extent is equal to this value and then it will be solved. This ensures the interaction potential has the same shape.
    end

    properties (SetAccess = private)
        potentialParameters
        ScaleInvarient_Potential_Minimum_Location
    end

    properties (SetAccess = private, Hidden)
        potentialParameterIdx = [1 1 1];
        InteractionOptions = struct('type','SRALRR','params',[]);
        
        dist = 'euc';
        dist_arg = [];
        
        Data_Range = []; % Dependent property set when Data_Storage_Centers is set. This gives the range of the data in each dimension.
    end

    properties (Hidden)
        Scale_Factor = 1; % This changes for each object and is not a global option
        Object_Scale = 1; % This changes for each object and is not a global option

        Use_ConvexHull = true;
        
        % These next two parameters would make apearence in functions
        % computing boundary curvature.
        Curvature_Smoothing_Size     = 2; % Not used
        Curvature_Max_Radius         = 35; % Not used 
        
        Potential_PreMultiplier = 1; % Not used
    end

    methods
        function options = seedPointOptions(varargin)

            % First get the potential paramters
            try
                options.potentialParameters = load('potentialParameters.mat');
            catch
                error('seedPointOptions:missingPotentialParamters','The file potentialParameters.mat is missing or not on the path. This file is required for setting the potentialParamters. It may be created with the function computePotentialParamters if you do not have it.')
            end

            depth_idx = 1;
            center_idx = 1;
            extent_idx = numel(options.potentialParameters.extent);

            options.potentialParameterIdx = [depth_idx, center_idx, extent_idx];
            idx = options.potentialParameters.parametersIdx == encodePotentialIdx(depth_idx, center_idx, extent_idx);
            options.InteractionOptions.params = options.potentialParameters.parameters(idx,:);

            options.Potential_Depth = options.potentialParameters.depth(1);
            options.Potential_Minimum_Location = options.potentialParameters.center(1);
            options.Potential_Extent = options.potentialParameters.extent(end);

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

        function obj = set.Solver_Space_Scale(obj, value)
            if ischar(value)
                if ~strcmp(value,'auto')
                    error('seedPointOptions:badSolverScale','Solver_Space_Scale must be set to ''auto'' or a positive finite real scalar value.')
                end
            else
                validateattributes(value,{'numeric'},{'positive','scalar','real','finite'})
                value = double(value);
                
                if value < 10
                    warning('seedPointOptions:potentiallyBadSolverScale','The value of Solver_Space_Scale, %f, may be small and lead to problems with the interaction potential',value)
                end
            end
            
            obj.Solver_Space_Scale = value;            
        end
        
        function obj = set.Data_Storage_Centers(obj, value)
            validateattributes(value,{'cell'})
            obj.Data_Storage_Centers = value;
            
            dat_range = cellfun(@(x) range(x), value);
            obj.Data_Range = dat_range(:)'; %#ok<MCSUP>
        end
        
        function obj = set.Wigner_Seitz_Radius(obj,value)
            validateattributes(value,{'double'},{'positive','scalar','real','finite'})
            obj.Wigner_Seitz_Radius = value;
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

        function obj = set.Potential_Depth(obj,value)
            validateattributes(-value,{'double'},{'scalar','nonnegative','real','finite'})

            if ~isempty(obj.Potential_Minimum_Location) && ~isempty(obj.Potential_Extent) %#ok<MCSUP>
                obj = setPotentialParameter(obj, value, obj.Potential_Minimum_Location, obj.Potential_Extent); %#ok<MCSUP>
            end
            obj.Potential_Depth = value;
        end

        function obj = set.Potential_Minimum_Location(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})

            if value > obj.Potential_Extent %#ok<MCSUP>
                error('seedPointOptions:badSet','Potential_Minimum_Location must be smaller than Potential_Extent.')
            end

            if ~isempty(obj.Potential_Depth) && ~isempty(obj.Potential_Extent) %#ok<MCSUP>
                obj = setPotentialParameter(obj, obj.Potential_Depth, value, obj.Potential_Extent); %#ok<MCSUP>
            end
            obj.Potential_Minimum_Location = value;
        end

        function obj = set.Potential_Extent(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})

            if value < obj.Potential_Minimum_Location %#ok<MCSUP>
                error('seedPointOptions:badSet','Potential_Extent must be larger than Potential_Minimum_Location.')
            end
            if ~isempty(obj.Potential_Depth) && ~isempty(obj.Potential_Minimum_Location) %#ok<MCSUP>
                obj = setPotentialParameter(obj, obj.Potential_Depth, obj.Potential_Minimum_Location, value); %#ok<MCSUP>
            end
            obj.Potential_Extent = value;
        end

        function obj = set.Potential_Padding_Size(obj,value)
            validateattributes(value,{'double'},{'integer','scalar','nonnegative','real','finite'})
            obj.Potential_Padding_Size = value;
        end

        function obj = set.Maximum_Initial_Potential(obj,value)
            validateattributes(value,{'double'},{'scalar','real','>',-Inf,'<=',Inf})
            if value < obj.Minimum_Initial_Potential %#ok<MCSUP>
                error('seedPointOptions:badInput','MaximumM_Initial_Potential must be larger than inimum_Initial_Potential.')
            end
            obj.Maximum_Initial_Potential = value;
        end
        
        function obj = set.Minimum_Initial_Potential(obj,value)
            validateattributes(value,{'double'},{'scalar','real','>',-Inf,'<=',Inf})
            if value > obj.Maximum_Initial_Potential %#ok<MCSUP>
                error('seedPointOptions:badInput','Minimum_Initial_Potential must be smaller than Maximum_Initial_Potential.')
            end
            obj.Minimum_Initial_Potential = value;
        end

        function obj = set.Potential_Scale(obj,value)
            validateattributes(value,{'double'},{'scalar','real','positive','nonzero'})
            if isinf(value)
                error('setPotentialScale:inf','Input must not be Inf.')
            end
            obj.Potential_Scale = value;
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

        function obj = set.Minimum_Hole_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','real','finite'})
            obj.Minimum_Hole_Size = value;
        end

        function obj = set.Curvature_Smoothing_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','integer','real','finite'})
            obj.Curvature_Smoothing_Size = value;
        end

        function obj = set.Curvature_Max_Radius(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            obj.Curvature_Max_Radius = value;
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
            obj.Use_Parallel = value;
        end

        function obj = set.Maximum_Memory(obj,value)
            validateattributes(value,{'double'},{'positive','scalar','real'})
            obj.Maximum_Memory = value;
        end
        
        function obj = set.Debug(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('seedPointOptions:badInput','Expected input to be logical.')
            end
            obj.Debug = value;
        end

        function obj = set.Object_Of_Interest(obj,value)
            if isempty(value)
                obj.Object_Of_Interest = value;
            else
                validateattributes(value,{'double'},{'scalar','integer','nonnegative','real','finite'})
                obj.Object_Of_Interest = value;
            end
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

        function obj = set.Scale_Factor(obj,value)
            validateattributes(value,{'double'},{'scalar','real','finite','positive'})
            obj.Scale_Factor = value;
        end

        function obj = set.ScaleInvarient_Potential_Extent(obj, value)
            validateattributes(value,{'double'},{'scalar','real','finite','positive'})
            obj.ScaleInvarient_Potential_Extent = value;
            obj = setPotentialParameter(obj,obj.Potential_Depth,obj.Potential_Minimum_Location,obj.Potential_Extent); %#ok<MCSUP>
        end
        
        function plotPotential(obj)

            % Compute the interaction potential
            scaleFactor = obj.Potential_Extent / obj.ScaleInvarient_Potential_Extent;
            
            x = obj.InteractionOptions.params;
            r = 0:0.05:1.3*obj.ScaleInvarient_Potential_Extent;
            Vint = 1./(r+0.2) - x(1)*exp(-(r-x(2)).^2/(2*x(3)^2));
            r = r * scaleFactor;
            
            figure
            
            % Plot interaction potential
            ax1 = subplot(2,1,1);
            line([r(1),r(end)], [0 0],'linestyle','--','color','k')
            line(r,Vint,'color','b','linewidth',2)
            
            xticks = [0, obj.Potential_Minimum_Location, obj.Potential_Extent];
            set(gca,'XTick',xticks,'XLim',[r(1),r(end)],'YLim',1.2*abs(obj.Potential_Depth)*[-1,1]);
            title(sprintf('Interaction potential\n (d_0=%0.2f, r_0=%0.2f, r_a=%0.2f) @ r_{a,SI}=%0.2f', obj.Potential_Depth, obj.Potential_Minimum_Location, obj.Potential_Extent, obj.ScaleInvarient_Potential_Extent))
            
            % Plot interaction force
            ax2 = subplot(2,1,2);
            line([r(1),r(end)], [0 0],'linestyle','--','color','k')
            line(r(1:end-1)+0.025,diff(Vint),'color','b','linewidth',2)
            
            xticks = [0, obj.Potential_Minimum_Location, obj.Potential_Extent];
            set(gca,'XTick',xticks,'XLim',[r(1),r(end)],'YLim',1.2*max(diff(Vint))*[-1,1]);
            title('Interaction force')
            
            % Set theme
            setTheme(gcf,'light')

        end
    end

    methods (Access = private)
        function obj = setPotentialParameter(obj,depth,center,extent)

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
            % TODO: Add function that checks to see that the attractive
            % extent, depth, and center actually are where we set them to
            % be. Not all combinations of parameters are possible (small
            % center and large extent for example)
        end
        
        function validateInteractionPotential(obj)
            % Confirm that the interaction potential has the correct depth,
            % minimum location, and extent; give a warning if the values
            % are wrong.
            %
            % Wrong values occur when A, mu, and sigma could not be solved
            % for to make the the depth, minimum location, and extent the
            % same as the values set.
            
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
            y = obj.Potential_Depth;
            y_hat = Vint(y_hat);
            y_err = abs(y_hat - y)/y;
            
            if y_err > 1e-2
                warning('seedPointOptions:unsolvableInteractionPotential', 'The potential depth is %0.2f%% different than the set potential depth. Consider modifying the potential parameters to ensure the interaction potential is as expected.',y_err*100)
            end
        end
    end
end

function out = encodePotentialIdx(depth_idx,center_idx,extent_idx)
    out = uint32(depth_idx + bitshift(center_idx,8) + bitshift(extent_idx,16));
end

