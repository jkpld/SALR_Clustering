function [centers, varargout] = computeObjectCenters(BW,B,intraS,initialLocations,options)
% COMPUTEOBJECTCENTERS  Compute object centers by particle simulation
%
% centers = computeObjectCenters(BW,B,markers,options)
% [centers, Info] = computeObjectCenters(BW,B,markers,options)
%
% Input parameters:
%
% BW - Object mask
%
% B - Object boundary. Multiple contours should be delimeted from each
%     other by row of NaN's. Outer contours should be clockwise oriented,
%     inner contours (holes) should be counter-clockwise oriented.
%
% initialLocations - The set of possible initial locations. Each row is a
%                    locations, each column is a dimension (x,y).
%
% options - element of class declumpOptions
%
% Output parameters:
%
% centers - Locations of the centers found. Each row is a location, each
%           column a dimension.
%
% Info - If options.Debug is set to true, then Info will be a structure
%        array with the following fields.
%
%        r0 - Initial particle locations (rows are locations)
%        ComputeIntitialPoints - Information structure returned by
%                                computeInitialPoints function.
%        V - Confinement potential
%        dVx - x derivative of V
%        dVy - y derivative of V
%        r_end - Final particle locations after simulation
%        solverTime - Computation time of solution
%        ode_solution - ode solution structure returned by ode23
%        converged - logical. true if the solution converged before
%                    stopping at maximum time
%
% See also COMPUTEINITIALPOINTS INTERACTINGPARTICLESYSTEM
% INTERACTINGPARTICLESYSTEM_CONVERGEEVENT

% James Kapaldo

WIGNER_SEITZ_RADIUS     = options.Wigner_Seitz_Radius;
INITIAL_SPEED           = options.Initial_Speed;
InteractionOptions      = options.InteractionOptions;

POINT_SELECTION_METHOD  = options.Point_Selection_Method;

PARTICLE_DAMPING_RATE   = options.Particle_Damping_Rate;
SOLVER_TIME_RANGE       = options.Solver_Time_Range;

CHARGE_NORMALIZATION_BETA = options.Charge_Normalization_Beta;

PAD_SIZE                = 5;

debug = options.Debug;

% Process the mask -------------------------------------------------------

% Get mask size
maskSize = size(BW);

% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*[1 1]);
intraS = padarray(intraS,PAD_SIZE*[1 1],1);



% Compute confining potential --------------------------------------------

% Get potential inside well. Compute the binary distance transform, smooth
% it out, and invert it to get 1/r.
V = double(bwdist(~BW_pad));
V = imfilter(V,fspecial('gaussian',7,1));
V = 1./V;

% Set potential outside well to be strongly increasing with distance. Don't
% set to infinity because, when modeling, if a particle happened to go into
% the Inf region during a time step, then the solution will break.
V_out = (bwdist(BW_pad)+1);
V_out = V_out.^(V_out+1);
V(~BW_pad) = V_out(~BW_pad);

Vdist = V;
V = V.*intraS;

% Compute gradient of confining potential --------------------------------
[dVx,dVy] = gradient(V); 

if debug
    Info.Vdist = Vdist(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
    Info.V = V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
    Info.dVx = dVx(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
    Info.dVy = dVy(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);
end


% Compute particle initial positions -------------------------------------

% This is computed near the top of this function because if there is only 1
% particle, then we will not run declumping on this object.

% Particle locations
Options.rs = WIGNER_SEITZ_RADIUS;
Options.curvatureCenters = initialLocations;



if debug
    [r0,ComputInitalPointsInfo] = computeInitialPoints(POINT_SELECTION_METHOD,BW,B,V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE),Options);
    Info.ComputeInitialPoints = ComputInitalPointsInfo;
    Info.r0 = r0;
else
    r0 = computeInitialPoints(POINT_SELECTION_METHOD,BW,B,V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE),Options);
end

% Number of particles
N = size(r0,1);

% If the number of particles is 1, then do continue further to model.
% Return with zero centers, this object does not need to be declumped.
if N < 2 || ~isfinite(N)
    centers = [];
    if debug
        Info.Vdist = [];
        Info.V = [];
        Info.dVx = [];
        Info.dVy = [];
        Info.r_end = [];
        Info.solverTime = nan;
        Info.ode_solution = [];
        Info.converged = 1;
        varargout{1} = Info;
    end
    return;
end



% Create interpolating functions for confining force and potential -------

% This method is slower --- (see below for faster interpolation)
% [Y,X] = ndgrid(1:size(BW_pad,1),1:size(BW_pad,2));
% V = griddedInterpolant(Y,X,V);
% dVx = griddedInterpolant(Y,X,dVx);
% dVy = griddedInterpolant(Y,X,dVy);
% Note that the above interpolation methods will use bilinear interpolation
% with bilinear extrapolation by default.

% Use faster interpolation methods with less overhead. 
dVx = @(r) interp2mex(dVx,r(:,2),r(:,1)); 
dVy = @(r) interp2mex(dVy,r(:,2),r(:,1));

% These mex interpolation functions assume the input is from an image (x
% and y spacing of 1) and it uses nearest neighbor extrapolation.

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

% Initialize particles ---------------------------------------------------

% Add in the padding to the initial locations.
r0 = r0 + PAD_SIZE;

% Particle initial velocities.
v0 = rand(N,2);
v0 = INITIAL_SPEED * v0 ./ sqrt(sum(v0.^2,2));

% Particle charge and mass
q = N^(-CHARGE_NORMALIZATION_BETA);
m = 1;

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

% Set up inputs to interactingParticleSystem()
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

% Create the funciton to solve
odeFun = @(t,y) interactingParticleSystem(t,y,SystemInputs);

% Set up the solver options.
ode_options = odeset('Events',@interactingParticleSystem_convergeEvent,...
    'Vectorized','on',...
    'RelTol',1e-4,...
    'AbsTol',1e-6,...
    'NormControl','on');

% Solve the system -------------------------------------------------------

% Record the time it takes to solve.
tic

sol = ode23(odeFun,SOLVER_TIME_RANGE,y0,ode_options);

solverTime = toc;

if isempty(sol.xe) || sol.xe == SOLVER_TIME_RANGE(end)
    Info.converged = 0;
    warning('ComputeObjectCenters:noConverge','The ODE solver did not converge in the time range given.')
else
    Info.converged = 1;
end

% Extract final solution and find centers of clusters --------------------

% Final solution
y_end = sol.y(:,end);

% Particle pixel locations.
r = round(y_end(rInds)) - PAD_SIZE; % Subtract the PAD_SIZE

% It could be, if the solution did not converge, that there are particle
% positions outside of the image. Remove all of these values.
r(r(:,1) < 1 | r(:,1) > maskSize(1) | r(:,2) < 1 | r(:,2) > maskSize(2),:) = [];

% Linear indices of particle locations
rInd = r(:,1) + (r(:,2)-1) * maskSize(1);

rInd(~isfinite(rInd)) = [];

% Search for the clusters by creating an image with the particles as 1's,
% then dilate and find the centroids of the connected components.

mask = false(maskSize);
mask(rInd) = 1;
mask = imdilate(mask,strel('disk',ceil(options.Potential_Minimum_Location)));

CC = bwconncomp(mask);
props = regionprops(CC,'Centroid'); % Perhaps should replace this with faster code for getting centroids.

% Combine all of the centroids into one array and flip the xy components as
% regionprops changes the order.
centers = fliplr(cat(1,props.Centroid));

% Add the top left corner of the boundary back in
% centers = centers + minB; 

if debug
    Info.r_end = r;
    sol.y(rInds,:) = sol.y(rInds,:) - PAD_SIZE; % Offset all solution positions by PAD_SIZE, to help preventing confusion later on.
    Info.solverTime = solverTime;
    Info.ode_solution = sol;
    varargout{1} = Info;
end

end