classdef displayProgress < handle
    % DISPLAYPROGRESS Class to help with displaying progress during
    % calculation. 
    %
    % progress = displayProgress(num_iterations, num_updates, active, is_parallel)
    % 
    % If Matlab version is R2017a or newer, then the progress can be
    % displayed even when using parallel computing. If the Matlab version
    % is not this new, then progress cannot be monitored when using
    % parallel computing.
    properties (SetAccess = public)
        active = true
        is_parallel = false
        number_of_iterations
        name
    end

    properties (SetAccess = private)
        iteration
        timer
        generate_display_at
    end
    
    properties (SetAccess = private, Hidden)
        first_print = 0;
        num_chars_written = 0;
    end

    methods
        function obj = displayProgress(varargin)
            % obj = displayProgress(number_of_iterations, number_of_displays, active)

            p = inputParser;
            p.FunctionName = 'displayProgress';

            addRequired(p,'number_of_iterations', @(t) validateattributes(t,{'numeric'},{'integer','scalar','positive'}))
            addOptional(p,'number_of_displays', 5, @(t) validateattributes(t,{'numeric'},{'integer','scalar','positive'}))
            addOptional(p,'active', obj.active, @(t) t==0 || t==1)
            addOptional(p,'is_parallel', obj.is_parallel, @(t) t==0 || t==1)
            addParameter(p,'name','Iteration',@(t) ischar(t))
            parse(p,varargin{:})

            N = p.Results.number_of_iterations;
            K = p.Results.number_of_displays;
            active = p.Results.active;
            is_parallel = p.Results.is_parallel;
            obj.name = p.Results.name;

            obj.number_of_iterations = N;
            obj.iteration = 0;
            
            showAt = unique(round(linspace(1,N,K+1)));
            if numel(showAt) > 1
                showAt(1) = [];
            end
            
            obj.generate_display_at = showAt;
            obj.active = logical(active);
            obj.is_parallel = is_parallel;
        end
        
        function set.is_parallel(obj,value)
            if value~=0 && value~=1
                error('displayProgress:invalidValue','is_parallel should be 0 or 1 or logical.')
            end
            obj.is_parallel = value;
            
            if obj.active && value && verLessThan('matlab','9.2') %#ok<MCSUP>
                % Verbose output with parallel computing not available with
                % matlab versions less than 9.2.
                warning('displayProgress:parallelOutputNotAvailable','Verbose output with parallel computing not available with matlab versions less than 9.2.')
                obj.active = false; %#ok<MCSUP>
            end
        end
        
        function set.active(obj, value)
            validateattributes(value,{'logical'},{'scalar'})
            obj.active = value;
            
            if obj.is_parallel && obj.active && verLessThan('matlab','9.2') %#ok<MCSUP>
                % Verbose output with parallel computing not available with
                % matlab versions less than 9.2.
                warning('displayProgress:parallelOutputNotAvailable','Verbose output with parallel computing not available with matlab versions less than 9.2.')
                obj.active = false;
            end
        end
        
        function Que = start(obj)
            Que = [];
            if ~obj.active
                return;
            end
            if obj.is_parallel
                Que = parallel.pool.DataQueue;
                afterEach(Que, @obj.iteration_end);
            end
            obj.timer = tic();
        end

        function iteration_end(obj,varargin)
            if ~obj.active
                return;
            end
            
            obj.iteration = obj.iteration + 1;
            
            if obj.iteration > obj.number_of_iterations
                obj.active = false;
                obj.timer = [];
                return;
            end
            
            currentTime = toc(obj.timer)/60;
            expectedTime = currentTime * obj.number_of_iterations / obj.iteration;

            if any(obj.iteration == obj.generate_display_at)
                str = sprintf('%s >> %s %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),obj.name,obj.iteration,obj.number_of_iterations,currentTime,expectedTime);
%                 fprintf('%s', str);
                
                % Comment out the above fprintf statement and uncomment the
                % below code to have iteration lines print on the same
                % line. (At a new iteration we print out enough backspaces
                % to remove the previously printed line, then we print the
                % new line.) The code works well as long as there is no
                % other output (like error messages) during the iterations.
                
                if obj.first_print == 0
                    obj.num_chars_written = length(str);
                    obj.first_print = 1;
                    fprintf('%s', str);
                else
                    fprintf(repmat('\b',1,obj.num_chars_written))
                    fprintf('%s', str);
                    obj.num_chars_written = length(str);
                end
            end
        end
    end
end


%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
