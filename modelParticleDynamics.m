function [r, varargout] = modelParticleDynamics(V,r0,options)
% MODELPARTICLEDYNAMICS Model particles confined in potential V with
% initial positions r0 and return their approximate equilibrium positions.
%
% r = modelParticleDynamics(V,r0,options)
% [r, Info] = modelParticleDynamics(V,r0,options)
%
% Input parameters:
%
% V         -   Confining potential
% r0        -   N x d array giving the initial position of the N particles
%                 in d-dimensions
% options   -   Element of class seedPointOptions
%
% Output parameters:
%
% r    -   Final particle positions
% Info -   If options.Debug is true, then Info will be a structure
%          with the following fields.
%
%             dVx - x derivative of V
%             dVy - y derivative of V
%             ode_solution - ode solution strucutre returned by
%                            ode23
%             converged - logical flag, if true, then the solution
%                         converged before stopping at maximum time
%
% See also EXTRACTCLUSTERCENTERS COMPUTEOBJECTSEEDPOINTS PROCESSOBJECTS

% James Kapaldo

INITIAL_SPEED               = options.Initial_Speed;
PARTICLE_DAMPING_RATE       = options.Particle_Damping_Rate;
SOLVER_TIME_RANGE           = options.Solver_Time_Range;
CHARGE_NORMALIZATION_BETA   = options.Charge_Normalization_Beta;
PAD_SIZE                    = options.Potential_Padding_Size;
DEBUG                       = options.Debug;
InteractionOptions          = options.InteractionOptions;
MASS_CHARGE_MULTIPLIER      = options.Mass_Charge_Multiplier;

MAX_FAILS = 5;
RESTART_UNDER = 5;
ON_FAIL_BACKTRACK = 2;
HISTORY_SIZE = 5;

% Number of particles
N = size(r0,1);

% Object image size
maskSize = size(V) - 2*PAD_SIZE;

% Compute gradient of confining potential --------------------------------
[dVx,dVy] = gradient(V);

% if DEBUG
%     Info.dVx = dVx(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
%     Info.dVy = dVy(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
% end

% Create interpolating functions for confining force and potential -------
dVx = @(r) interp2mex(dVx,r(:,2),r(:,1));
dVy = @(r) interp2mex(dVy,r(:,2),r(:,1));

% These mex interpolation functions assume the input is from an image (x
% and y spacing of 1) and it uses nearest neighbor extrapolation.

% .......................................................................
% Note, if the above mex functions are not installed, then the code here
% can be used (though it is slower)
%
% [Y,X] = ndgrid(1:size(BW_pad,1),1:size(BW_pad,2));
% V = griddedInterpolant(Y,X,V);
% dVx = griddedInterpolant(Y,X,dVx);
% dVy = griddedInterpolant(Y,X,dVy);
% Note that the above interpolation methods will use bilinear interpolation
% with bilinear extrapolation by default.
% .......................................................................


% Create particle interaction potential and force ------------------------
switch InteractionOptions.type
    case 'Coulomb'
        % Vint  = @(D)  1./D;    % Coulomb potential
        dVint = @(D) -1./D.^2; % Coulomb force
    case 'SRALRR'
        % 1/r repulsion with gaussian attraction
        A = InteractionOptions.params(1);
        mu = InteractionOptions.params(2);
        sig = InteractionOptions.params(3);

        % Vint  = @(D)  1./D - A * exp(-(D-mu).^2/(2*sig^2));
        dVint = @(D) -1./(D+0.2).^2 + (A*(D-mu)/(sig^2)) .* exp(-(D-mu).^2/(2*sig^2));
end

% Add in the padding to the initial locations.
r0 = r0 + PAD_SIZE;

% Particle initial velocities.
v0 = rand(N,2);
v0 = INITIAL_SPEED * v0 ./ sqrt(sum(v0.^2,2));

% Particle charge and mass
q = sqrt(MASS_CHARGE_MULTIPLIER) * N^(-CHARGE_NORMALIZATION_BETA);
m = 1/MASS_CHARGE_MULTIPLIER;

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
% create an array that is N*(N-1)/2 x 2. These unit vectors are simply
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
% SystemInputs.d = size(r0,2);
SystemInputs.q = q;
SystemInputs.m = m;
SystemInputs.alpha = alpha;
SystemInputs.dVx = dVy; % Note that the order of dVx dVy is switched. This is because y is the first dimension and x is the second.
SystemInputs.dVy = dVx;
SystemInputs.dVint = dVint;
SystemInputs.N_to_NtNm1o2 = N_to_NtNm1o2;
SystemInputs.pdistInds = pdistInds;

% Put initial conditions in correct form
offset = (0:4:N*4-1)';

rInds = [1,2] + offset;
vInds = [3,4] + offset;

y0 = zeros(N,1);
y0(rInds) = r0;
y0(vInds) = v0;

% Create the function to solve
odeFun = @(t,y) interactingParticleSystem(t,y,SystemInputs);

% Initialize ode history
hstry = ode_history(SOLVER_TIME_RANGE(1), y0, HISTORY_SIZE);

% Set up the solver options.
ode_options = odeset('Events',@(t,y) interactingParticleSystem_convergeEvent(t,y,m,hstry),...
    'Vectorized','on',...
    'RelTol',1e-4,...
    'AbsTol',1e-6,...
    'NormControl','on');%,...
%     'InitialStep',InitialStep);


% Solve the system -------------------------------------------------------
% Record the time it takes to solve.
startClock = tic;

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

        failIdx = find(any(isnan(sol.y),1),1);

        if failIdx < RESTART_UNDER
%             fprintf('start fail %d\n',quitIterations+1)
            ode_options.InitialStep = (sol.x(2)-sol.x(1))/2;
            hstry.hardreset(SOLVER_TIME_RANGE(1), y0);
        else
%             fprintf('large jump fail\n')
            returnTo = failIdx - ON_FAIL_BACKTRACK;
            startTime = sol.x(returnTo);

            ode_options.InitialStep = (sol.x(returnTo+1) - startTime)/2;

            sol_hist.y = [sol_hist.y, sol.y(:,1:returnTo-1)];
            sol_hist.x = [sol_hist.x, sol.x(1:returnTo-1)];

            hstry.rewrite(sol.x(returnTo-(HISTORY_SIZE):returnTo-1), sol.y(:,returnTo-(HISTORY_SIZE):returnTo-1), ON_FAIL_BACKTRACK)

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

% Particle pixel locations.
r = round(y_end(rInds)) - PAD_SIZE; % Subtract the PAD_SIZE

% It could be, if the solution did not converge, that there are particle
% positions outside of the image. Remove all of these values.
r(r(:,1) < 1 | r(:,1) > maskSize(1) | r(:,2) < 1 | r(:,2) > maskSize(2),:) = [];


if DEBUG
    Info.solverTime = solverTime;
    Info.ode_solution = sol;
    Info.SystemInputs = SystemInputs;
    varargout{1} = Info;
end

end
