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
%   (1, Info) , [pixels]
%
% Maximum_Initial_Potential - The maximum confining potential value that a
%   particle can have at its initial position. Any possible initial
%   position with a confining potential larger than this value will not be
%   used as an initial particle position. 
%   (-Inf, Inf] , [arb. units]
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
% Curvature_Smoothing_Size - The standard deviation of the gaussians used
%   for smoothing the boundary and calculating the curvature.
%   [0, Inf) & integer , [pixels]
%
% Use_GPU - Determines if a GPU will be used to speed up calculation.
%   logical
%
% Use_Parallel - Determines if multiple CPUs will be sued to speed up
%   calculation.
%   logical
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
    
    properties

        Wigner_Seitz_Radius          = 10;
        Potential_Depth
        Potential_Minimum_Location
        Potential_Extent
        Potential_Padding_Size       = 5;
        Maximum_Initial_Potential    = 1/3;
        Initial_Speed                = 0.01;
        Particle_Damping_Rate        = 5e-4;
        Charge_Normalization_Beta    = 1/3;
        Solver_Time_Range            = 0:10:1500;
        Point_Selection_Method       = 'r0set_uniformRandom';
        
        Minimum_Hole_Size            = 50;
        
        Curvature_Smoothing_Size     = 2;
        
        Use_GPU                      = false;
        Use_Parallel                 = false;
        
        Debug                        = false;
        
        Object_Of_Interest           = [];
    end
    
    properties (SetAccess = private)
        potentialParameters
    end
    
    properties (SetAccess = private, Hidden)
        potentialParameterIdx = [1 1 1];
        InteractionOptions = struct('type','SRALRR','params',[]);
    end
    
    methods
        function options = seedPointOptions(varargin)
            
            % First get the potential paramters
            try
                options.potentialParameters = load('potentialParameters.mat');
            catch
                error('seedPointOptions:missingPotentialParamters','The file potentialParameters.mat is missing or not on the path. This file is required for setting the potentialParamters. It may be created with the function computePotentialParamters if you do not have it.')
            end
            
            options.Potential_Depth = options.potentialParameters.depth(1);
            options.Potential_Minimum_Location = options.potentialParameters.center(1);
            options.Potential_Extent = options.potentialParameters.extent(end);
            
            options.InteractionOptions.params = permute(options.potentialParameters.parameters(1,1,end,:),[4,1,2,3]);
            
            
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
            
            idx = find(obj.potentialParameters.depth == value); %#ok<MCSUP>
            
            if isempty(idx)
                badValue = value;
                [~, idx] = min(abs(obj.potentialParameters.depth-value));  %#ok<MCSUP>
                value = obj.potentialParameters.depth(idx);  %#ok<MCSUP>
                warning('seedPointOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Depth %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
            end
            
            obj.potentialParameterIdx(1) = idx; %#ok<MCSUP>
            obj.InteractionOptions.params = permute(obj.potentialParameters.parameters(obj.potentialParameterIdx(1),obj.potentialParameterIdx(2),obj.potentialParameterIdx(3),:),[4,1,2,3]); %#ok<MCSUP>
            obj.Potential_Depth = value;
        end
        
        function obj = set.Potential_Minimum_Location(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            
            idx = find(obj.potentialParameters.center == value); %#ok<MCSUP>
            
            if isempty(idx)
                badValue = value;
                [~, idx] = min(abs(obj.potentialParameters.center-value));  %#ok<MCSUP>
                value = obj.potentialParameters.center(idx);  %#ok<MCSUP>
                warning('seedPointOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Minimum_Location %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
            end
            
            if value > obj.Potential_Extent %#ok<MCSUP>
                error('seedPointOptions:badSet','Potential_Minimum_Location must be smaller than Potential_Extent.')
            end
            
            obj.potentialParameterIdx(2) = idx; %#ok<MCSUP>
            obj.InteractionOptions.params = permute(obj.potentialParameters.parameters(obj.potentialParameterIdx(1),obj.potentialParameterIdx(2),obj.potentialParameterIdx(3),:),[4,1,2,3]); %#ok<MCSUP>
            obj.Potential_Minimum_Location = value;
        end
        
        function obj = set.Potential_Extent(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            
            idx = find(obj.potentialParameters.extent == value); %#ok<MCSUP>
            
            if isempty(idx)
                badValue = value;
                [~, idx] = min(abs(obj.potentialParameters.extent-value));  %#ok<MCSUP>
                value = obj.potentialParameters.extent(idx);  %#ok<MCSUP>
                warning('seedPointOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Extent %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
            end
            
            if value < obj.Potential_Minimum_Location %#ok<MCSUP>
                error('seedPointOptions:badSet','Potential_Extent must be larger than Potential_Minimum_Location.')
            end
            
            obj.potentialParameterIdx(3) = idx; %#ok<MCSUP>
            obj.InteractionOptions.params = permute(obj.potentialParameters.parameters(obj.potentialParameterIdx(1),obj.potentialParameterIdx(2),obj.potentialParameterIdx(3),:),[4,1,2,3]); %#ok<MCSUP>
            obj.Potential_Extent = value;
        end
        
        function obj = set.Potential_Padding_Size(obj,value)
            validateattributes(value,{'double'},{'integer','scalar','positive','nonzero','real','finite'})
            obj.Potential_Padding_Size = value;
        end
        
        function obj = set.Maximum_Initial_Potential(obj,value)
            validateattributes(value,{'double'},{'scalar','real','>',-Inf,'<=',Inf})
            obj.Maximum_Initial_Potential = value;
        end
        
        function obj = set.Minimum_Hole_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','real','finite'})
            obj.Minimum_Hole_Size = value;
        end
        
        function obj = set.Curvature_Smoothing_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','integer','real','finite'})
            obj.Curvature_Smoothing_Size = value;
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
        
    end
    
end





































