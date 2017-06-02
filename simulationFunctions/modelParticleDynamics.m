function [r, varargout] = modelParticleDynamics(dV,r0,options)
% MODELPARTICLEDYNAMICS Model particles confined in potential with gradient
% dV and initial positions r0 and return their approximate equilibrium
% positions.
%
% r = modelParticleDynamics(dV,r0,options)
% [r, Info] = modelParticleDynamics(dV,r0,options)
%
% Input parameters:
%
% dV : Function handle taking in particle positions and returning the
%   confining force (the gradient of the confining potential)
% r0 : N x D array giving the initial position of the N particles.
% options : Element of class seedPointOptions
%
% Output parameters:
%
% r : Final particle positions
% Info : If options.Debug is true, then Info will be a structure with the
%   following fields.
%
%   ode_solution : ode solution structure returned by ode23
%   converged : logical flag, if true, then the solution converged before
%       stopping at maximum time
%   solverTime : time to reach convergence
%   SystemInputs : extra inputs for the solver, see the code for more
%       information about this -- to actually output this, uncomment the
%       code near the bottom; this takes a lot of space.
%
%
% See also EXTRACTCLUSTERCENTERS COMPUTEOBJECTSEEDPOINTS PROCESSOBJECTS

% James Kapaldo

% Pre-amble --------------------------------------------------------------
INITIAL_SPEED               = options.Initial_Speed;
MASS                        = options.Mass;
COUPLING_CONSTANT           = options.Coupling_Constant;
PARTICLE_DAMPING_RATE       = options.Particle_Damping_Rate;
CHARGE_NORMALIZATION_BETA   = options.Charge_Normalization_Beta;
InteractionOptions          = options.InteractionOptions;
SOLVER_TIME_RANGE           = options.Solver_Time_Range;
DISTANCE_METRIC             = options.dist;
DISTANCE_METRIC_ARGUMENT    = options.dist_arg;
DEBUG                       = options.Debug;

% History parameters
MAX_FAILS = 5;          % Maximum number of fails before quitting
RESTART_UNDER = 5;      % Restart entire simulation if less than RESTART_UNDER time steps have occured before a fail.
ON_FAIL_BACKTRACK = 2;  % Go back ON_FAIL_BACKTRACK time steps after fail occures (and if more than RESTART_UNDER time steps have already occured).
HISTORY_SIZE = 5;       % The number of past time steps stored.

% Number of particles, dimension, and potential size
N = size(r0,1);
D = size(r0,2);

% Create particle interaction potential and force ------------------------
switch InteractionOptions.type
    case 'Coulomb'
        % Vint  = @(D)  1./D;    % Coulomb potential
        dVint = @(D) -1./(D + 0.2).^2; % Coulomb force
    case 'SRALRR'
        % 1/r repulsion with gaussian attraction
        A = InteractionOptions.params(1);
        mu = InteractionOptions.params(2);
        sig = InteractionOptions.params(3);

        % Vint  = @(D)  1./D - A * exp(-(D-mu).^2/(2*sig^2));
        dVint = @(D) -1./(D + 0.2).^2 + (A*(D-mu)/(sig^2)) .* exp(-(D-mu).^2/(2*sig^2));
end

% Initialized particles --------------------------------------------------

% Particle initial velocities.
v0 = rand(N,D);
v0 = INITIAL_SPEED * v0 ./ sqrt(sum(v0.^2,2));

% Particle charge and mass
q = sqrt(COUPLING_CONSTANT) * N^(-CHARGE_NORMALIZATION_BETA);
m = MASS;

% Particle damping
alpha = @(t) zeros(N,1) + t * PARTICLE_DAMPING_RATE;


% Set up solver functions ------------------------------------------------

% Get the particle pair indices.
pdistInds = getPdistInds(N);

% Create a sparse matrix that will take us between the seperate N*(N-1)/2
% particle pairs and the N particles.
N_pair = N*(N-1)/2;

S1 = sparse(1:N_pair,pdistInds(:,1),1,N_pair,N,N_pair);
S2 = sparse(1:N_pair,pdistInds(:,2),1,N_pair,N,N_pair);

N_to_NtNm1o2 = S1 - S2;

% Example: we have N particles with positions r0. We need to form the unit
% vector between each pair of particles. All of these unit vectors will
% create an array that is N*(N-1)/2 x D. These unit vectors are simply
% created by
%
% ur = N_to_NtNm1o2 * r0;
% ur = ur ./ sqrt(sum(ur.^2,2));
%
% Now we will have an array of size N*(N-1)/2 x 1 with the force between
% each pair of particles. We need to add up all of the forces on a particle
% to result in a Nx1 array of forces. This is simply done by
%
% Force = - N_to_NtNm1o2' * (ur .* forces * q1*q2);
%
% Credit for the idea of using a sparse matrix :
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/274779

% Set up inputs to interactingParticleSystem()
SystemInputs.D = D;
SystemInputs.q = q;
SystemInputs.m = m;
SystemInputs.alpha = alpha;
SystemInputs.dV = dV;
SystemInputs.dVint = dVint;
SystemInputs.N_to_NtNm1o2 = N_to_NtNm1o2;
SystemInputs.pdistInds = pdistInds;

SystemInputs.dist = DISTANCE_METRIC;
SystemInputs.dist_arg = DISTANCE_METRIC_ARGUMENT;

% Put initial conditions in correct form
offset = (0:2*D:N*2*D-1)';

rInds = (1:D) + offset;
vInds = rInds + D;

y0 = zeros(N,1);
y0(rInds) = r0;
y0(vInds) = v0;

% Create the function to solve
odeFun = @(t,y) interactingParticleSystem(t,y,SystemInputs);

% Initialize ode history
hstry = ode_history(SOLVER_TIME_RANGE(1), y0, HISTORY_SIZE);

% Set up the solver options.
ode_options = odeset('Events',@(t,y) interactingParticleSystem_convergeEvent(t,y,m,hstry,D),...
    'Vectorized','on',...
    'RelTol',1e-4,...
    'AbsTol',1e-6,...
    'NormControl','on');%,...
%     'InitialStep',InitialStep);


% Solve the system -------------------------------------------------------

startClock = tic; % Record the time it takes to solve.

quitIterations = 0;
timeRange = SOLVER_TIME_RANGE;

sol_hist.y = [];
sol_hist.x = [];

while 1

    sol = ode23(odeFun,timeRange,y0,ode_options);

    if sol.ie == 2
        % Error do to particle position or momentum being NaN. Usually
        % caused by taking a time step that is too big.
        if quitIterations > MAX_FAILS
            break;
        end

        failIdx = find(any(isnan(sol.y),1),1); % The time step at which the solution first failed.

        if failIdx < RESTART_UNDER
            % Start completely over with a smaller initial time step

%             fprintf('start fail %d\n',quitIterations+1)
            ode_options.InitialStep = (sol.x(2)-sol.x(1))/2;
            hstry.hardreset(SOLVER_TIME_RANGE(1), y0);
        else
            % Go back ON_FAIL_BACKTRACK time steps and start again with a
            % smaller time step.

%             fprintf('large jump fail\n')
            returnTo = failIdx - ON_FAIL_BACKTRACK;
            startTime = sol.x(returnTo);

            ode_options.InitialStep = (sol.x(returnTo+1) - startTime)/2;

            sol_hist.y = [sol_hist.y, sol.y(:,1:returnTo-1)];
            sol_hist.x = [sol_hist.x, sol.x(1:returnTo-1)];

            start = max(returnTo - HISTORY_SIZE,1);
            stop = max(start,returnTo - 1);

            hstry.rewrite(sol.x(start:stop), sol.y(:,start:stop))
            y0 = sol.y(:,returnTo);

            timeRange = [startTime, 1500];
        end

        quitIterations = quitIterations + 1;
    else
        break;
    end
end

if quitIterations > 0
    sol.y = [sol_hist.y, sol.y];
    sol.x = [sol_hist.x, sol.x];
end

solverTime = toc(startClock);

if quitIterations > MAX_FAILS
    Info.converged = 0;
    warning('ComputeObjectCenters:solverFailed','The ODE solver failed to converge MAX_FAILS times, %d.', MAX_FAILS)
elseif isempty(sol.xe) || sol.xe == SOLVER_TIME_RANGE(end)
    Info.converged = 0;
    warning('ComputeObjectCenters:noConverge','The ODE solver did not converge in the time range given.')
else
    Info.converged = 1;
end


% Extract final solution -------------------------------------------------

% Final solution
y_end = sol.y(:,end);

% Final particle locations.
r = y_end(rInds);

if DEBUG
    Info.solverTime = solverTime;
    Info.ode_solution = sol;
    % Info.SystemInputs = SystemInputs; % This takes up lots of space and does not really give much information that cannot be found or computed by other means.
    varargout{1} = Info;
else
    varargout{1} = [];
end

end
