classdef ode_history < handle
    % ODE_HISTORY A class for storing an ODE solution at previous time
    % steps while solving. Use the history with an ODE 'Events' function to
    % catch more complex behavior than can be determined with only the
    % current solution.
    %
    % history = ode_history(t0, y0, history_size)
    %
    % Properties:
    %
    % solHistory - rolling array for storing the history. The latest (most
    % current) solution is stored in solHistory(:,end).
    %
    % time - rolling array for storing the time of each of the past
    % solutions in solHistory
    %
    % windowSize - the number of solutions to store. numel(time) =
    % windowSize, and size(solHisotry,2) = windowSize
    %
    % iteration - The number increases each time a new solution is added.
    %
    % Methods:
    %
    % obj = ode_history(t, sol, w), create history object by initializing
    % with the initial time, initial time, and the windowSize.
    %
    % add(obj, t, sol) - add a new solution at time t
    %
    % hardreset(obj, t, sol) - reset the history by setting iteration = 1
    % and solHistory(:,end) = sol, time(end) = t;
    %
    % rewrite(obj, ts, sols, offset) - completely rewrite the history using
    % time = ts, solHistory = sols, and iteration = iteration - offset;

    % James Kapaldo

    properties
        solHistory      = [];
        time            = [];
        iterationNumber = 1;
        windowSize      = 5;
        validEntries    = 1;
    end

    methods
        function obj = ode_history(t, sol, w)
            if nargin > 0
                obj.windowSize = w;
                if w > 0
                    obj.time = zeros(1,w);
                    obj.time(end) = t;
                    obj.solHistory = zeros([size(sol,1),w], 'like', sol);
                    obj.solHistory(:,end) = sol;
                end
            end
        end

        function add(obj,t,sol)

            obj.solHistory(:,1) = sol;
            obj.time(1) = t;

            if obj.windowSize > 1
                obj.solHistory = circshift(obj.solHistory,-1,2);
                obj.time = circshift(obj.time,-1,2);
            end

            obj.iterationNumber = obj.iterationNumber + 1;
            obj.validEntries = min(obj.windowSize, obj.validEntries+1);
        end

        function hardreset(obj,t, sol)
            obj.solHistory(:,end) = sol;
            obj.time(end) = t;
            obj.validEntries = 1;
        end

        function rewrite(obj,ts,sols)
            n = length(ts);

            obj.solHistory(:,(obj.windowSize-n+1):obj.windowSize) = sols;
            obj.time((obj.windowSize-n+1):obj.windowSize) = ts;
            obj.validEntries = n;
        end
    end
end
