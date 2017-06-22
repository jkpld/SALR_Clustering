classdef declumpOptions
% DECLUMPOPTIONS  Set options needed for declumping.
%
% Input can be structure array or parameter value pairs. Options not set
% will be given default values. To see the default values, look at the
% output with no inputs:
%   defaultValues = setDeclumpOptions()
%
% declumpOptions Properties:
%
% Max_Radius - The maximum distance a boundry vertex can be from a center
%   to be assigned to it.
%   (0, Inf) , [pixels]
%
% Min_Angle - The minimum *dot product* between the line formed between a
%   boundry vertex and a center and the unit normal (pointing inward) at 
%   the boundary vertex. -- The angle could be found by acos(Min_Angle).
%   [-1, 1] , [unitless]
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
%   {'random','uniform','uniformRandom','curvatureRandom','curvatureUniformRandom'}
%
% Min_Interior_Angle - The minimum interior angle for the triangulation of
%   the centers.
%   [0, 180] , [degrees]
%
% Max_Interior_Angle - The maximum interior angle for the triangulation of
%   the centers.
%   [0, 180] & > Min_Interior_Angle , [degrees]
%
% Search_Radius - The radius of the region over which each cut vertex is
%   optimized over.
%   [0, Inf) & integer , [pixels]
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
        Max_Radius                   = 35;
        Min_Angle                    = 0.5;
        
        Wigner_Seitz_Radius          = 10;
        Potential_Depth
        Potential_Minimum_Location
        Potential_Extent
        Initial_Speed                = 0.01;
        Particle_Damping_Rate        = 5e-4;
        Charge_Normalization_Beta    = 1/3;
        Solver_Time_Range            = 0:10:1500;
        Point_Selection_Method       = 'curvatureUniformRandom';
        
        Min_Interior_Angle           = 20;
        Max_Interior_Angle           = 110;
        
        Search_Radius                = 7;
        
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
        function options = declumpOptions(varargin)
            
            % First get the potential paramters
            try
                options.potentialParameters = load('potentialParameters.mat');
            catch
                error('declumpOptions:missingPotentialParamters','The file potentialParameters.mat is missing or not on the path. This file is required for setting the potentialParamters. It may be created with the function computePotentialParamters if you do not have it.')
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
                        error('declumpOptions:badInput','Input must be structure with properties as fiels, or property/value list with even number of elements.')
                    end
                    for i = 1:2:numel(varargin)
                        options.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        function obj = set.Max_Radius(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            obj.Max_Radius = value;
        end
        
        function obj = set.Min_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',-1,'<=',1,'real'})
            obj.Min_Angle = value;
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
            value = validatestring(value,{'random','uniform','uniformRandom','curvatureRandom','curvatureUniformRandom'});
            obj.Point_Selection_Method = value;
        end
        
        function obj = set.Potential_Depth(obj,value)
            validateattributes(-value,{'double'},{'scalar','nonnegative','real','finite'})
            
            idx = find(obj.potentialParameters.depth == value); %#ok<MCSUP>
            
            if isempty(idx)
                badValue = value;
                [~, idx] = min(abs(obj.potentialParameters.depth-value));  %#ok<MCSUP>
                value = obj.potentialParameters.depth(idx);  %#ok<MCSUP>
                warning('declumpOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Depth %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
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
                warning('declumpOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Minimum_Location %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
            end
            
            if value > obj.Potential_Extent %#ok<MCSUP>
                error('declumpOptions:badSet','Potential_Minimum_Location must be smaller than Potential_Extent.')
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
                warning('declumpOptions:undefinedPotentialParameterValues', 'There are no potential parameters defined for the Extent %g. Using the closest defined value %g.\nUse the function computePotentialParameters to compute the potential paramters for additional values.',badValue,value)
            end
            
            if value < obj.Potential_Minimum_Location %#ok<MCSUP>
                error('declumpOptions:badSet','Potential_Extent must be larger than Potential_Minimum_Location.')
            end
            
            obj.potentialParameterIdx(3) = idx; %#ok<MCSUP>
            obj.InteractionOptions.params = permute(obj.potentialParameters.parameters(obj.potentialParameterIdx(1),obj.potentialParameterIdx(2),obj.potentialParameterIdx(3),:),[4,1,2,3]); %#ok<MCSUP>
            obj.Potential_Extent = value;
        end
        
        function obj = set.Min_Interior_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',0,'<=',180})
            if value > obj.Max_Interior_Angle %#ok<MCSUP>
                error('declumpOptions:badSet','Min_Interior_Angle must be smaller than Max_Interior_Angle.')
            end
            obj.Min_Interior_Angle = value;
        end
        
        function obj = set.Max_Interior_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',0,'<=',180})
            if value < obj.Min_Interior_Angle %#ok<MCSUP>
                error('declumpOptions:badSet','Max_Interior_Angle must be larger than Min_Interior_Angle.')
            end
            obj.Max_Interior_Angle = value;
        end
        
        function obj = set.Search_Radius(obj,value)
            validateattributes(value,{'double'},{'scalar','integer','nonnegative','real','finite'})
            obj.Search_Radius = value;
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
                error('declumpOptions:badInput','Expected input to be logical.')
            end
            obj.Use_GPU = value;
        end
        
        function obj = set.Use_Parallel(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('declumpOptions:badInput','Expected input to be logical.')
            end
            obj.Use_Parallel = value;
        end
        
        function obj = set.Debug(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('declumpOptions:badInput','Expected input to be logical.')
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





































