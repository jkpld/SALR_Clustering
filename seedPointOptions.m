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
% Wigner_Seitz_Radius - The effective size of each particle in the
%   simulation. This sets the density of the particles used. The
%   approximate number of particles used in a particular object will be
%   area(object)/(pi* r_s^2), where r_s is the Wigner_Seitz_Radius.
%   (0, Inf) , [pixels]
%
% Potential_Depth - The depth of the potential.
%   (-Inf,0] , [arb. units]
%
% Potential_Minimum_Location - The location of the potential minimum.
%   (0, Inf) , [pixels]
%
% Potential_Extent - The radius at which the potential goes from attractive
%   to repulsive.
%   (0, Inf) & > Potential_Minimum_Location , [pixels]
%
% Potential_Padding_Size - The size of the padding to apply to an object
%   mask. This is important so that the potential is defined some distance
%   away from the object when running the simulation.
%   (1, Inf) & integer, [pixels]
%
% Maximum_Initial_Potential - The maximum confining potential value that a
%   particle can have as its initial position. Any possible initial
%   position with a confining potential larger than this value will not be
%   used as an initial particle position.
%   (-Inf, Inf] , [arb. units]
%
% Minimum_Initial_Potential - The minimum confining potential value that a
%   particle can have as its initial position. Any possible initial
%   position with a confining potential smaller than this value will not be
%   used as an initial particle position.
%   (-Inf, Inf] , [arb. units]
%
% Potential_Scale - Value used to set the depth of the confining potential
%   well.
%   (0, Inf) , [arb. units]
%
% Distance_Metric - The metric used for measuring distance between
%   particles. If metric is Minkowski, then an extra agrument can be given
%   with the exponent.
%   {'euclidean'; 'cityblock'; 'chebychev'; {'minkowski', exponent}}
%
% Initial_Speed - The initial speed of the particles in the simulation.
%   (-Inf, Inf) , [pixels/time]
%
% Particle_Damping_Rate - The rate at which the damping of the particles
%   increases with simulation time.
%   [0, Inf) , [1/time^2]
%
% Charge_Normalization_Beta - The beta exponent for charge normalization
%   based on number of particles.
%   (-Inf, Inf), [unitless]
%
% Solver_Time_Range - Range over which the particles are simulated.
%   [0, Inf) , [time]
%
% Point_Selection_Method - The method used to initialize the particle
%   locations.
%   {'random','uniform','uniformRandom','r0set_random','r0set_uniformRandom'}
%
% Minimum_Hole_Size - The minimum hole size allowed in the mask.
%   [0, Inf) , [pixels^2]
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
%   not take up memory. (It is assumed that the confining potential is of
%   class double.)
%   [0, Inf], Gb
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

        Wigner_Seitz_Radius          = 10;
        
        Potential_Depth
        Potential_Minimum_Location
        Potential_Extent
        
        Potential_Padding_Size       = 5;
        Maximum_Initial_Potential    = Inf;
        Minimum_Initial_Potential    = -Inf;

        Potential_Scale              = NaN;
        
        Distance_Metric              = 'euclidean';

        Initial_Speed                = 0.01;
        Particle_Damping_Rate        = 5e-4;
        Charge_Normalization_Beta    = 1/3;
        Solver_Time_Range            = 0:10:1500;
        Point_Selection_Method       = 'r0set_uniformRandom';

        Minimum_Hole_Size            = 50;

        Use_GPU                      = false;
        Use_Parallel                 = false;

        Maximum_Memory               = 1;
        
        Debug                        = false;

        Object_Of_Interest           = [];
    end

    properties (SetAccess = private)
        potentialParameters
    end

    properties (SetAccess = private, Hidden)
        potentialParameterIdx = [1 1 1];
        InteractionOptions = struct('type','SRALRR','params',[]);
        
        dist = 'euc';
        dist_arg = [];
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
        Mass_Charge_Multiplier = 1; % Not used
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

        
        
        function obj = set.Wigner_Seitz_Radius(obj,value)
            validateattributes(value,{'double'},{'positive','scalar','real','finite'})
            obj.Wigner_Seitz_Radius = value;
        end

        function obj = set.Initial_Speed(obj,value)
            validateattributes(value,{'double'},{'scalar','real','finite'})
            obj.Initial_Speed = value;
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
        
        function obj = set.Mass_Charge_Multiplier(obj,value)
            validateattributes(value,{'double'},{'scalar','real','positive','nonzero','finite'})
            obj.Mass_Charge_Multiplier = value;
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
            validateattributes(value,{'double'},{'scalar','real','finite','positive','nonzero'})
            obj.Scale_Factor = value;
        end

        function plotPotential(obj)

            x = obj.InteractionOptions.params;
            r = 0:0.05:1.3*obj.Potential_Extent;
            Vint = 1./(r+0.2) - x(1)*exp(-(r-x(2)).^2/(2*x(3)^2));

            figure
            plot(r,Vint,'b','linewidth',2)
            title(sprintf('Interaction potential (%0.2f,%0.2f,%0.2f)', obj.Potential_Depth, obj.Potential_Minimum_Location, obj.Potential_Extent))
        end
    end

    methods (Access = private)
        function obj = setPotentialParameter(obj,depth,center,extent)

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

%             obj.Potential_Depth = depth;
%             obj.Potential_Minimum_Location = center;
%             obj.Potential_Extent = extent;
            obj.potentialParameterIdx = [depth_idx, center_idx, extent_idx];

            % TODO: Add function that checks to see that the attractive
            % extent, depth, and center actually are where we set them to
            % be. Not all combinations of parameters are possible (small
            % center and large extent for example)
        end
    end

end

function out = encodePotentialIdx(depth_idx,center_idx,extent_idx)
    out = uint32(depth_idx + bitshift(center_idx,8) + bitshift(extent_idx,16));
end
