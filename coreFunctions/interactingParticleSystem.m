function dy = interactingParticleSystem(t,y,extraInputs)
% INTERACTINGPARTICLESYSTEM Compute the change in particle states dy=[dr/dt,
% dp/dt] from the current states y=[r, p].
%
% dy = interactingParticleSystem(t,y,extraInputs)
%
% The parameter vector y should have the form
% y  = [ y1, y2, y3, ... , yn ]         % D*N x n matrix
% yi = [ r1; p1; r2; p2; ...; rN; pN ]  % D*N x 1 array
% ri = [ ri_x, ri_y, ..., ri_D ]'       % D x 1 array
% pi = [ pi_x, pi_y, ..., pi_D ]'       % D x 1 array
%
% Where N is the number of particles and where if n > 1 the ode solver is
% trying to compute multiple points at once.
%
% ExtraInputs : structure with the fields
%       D : number of dimensions
%       q : particle charges
%       m : particle masses
%       alpha : particle damping coefficients
%       dV : confining potential gradient
%       dVind : derivative of particle interaction potential
%       pdistInds : indices of particle pairs
%       N_to_NtNm1o2 : ( read as, N to N*(N-1)/2 ) sparse matrix that can
%           create particle pair vectors from particle positions and
%           accumulate all of the forces on a particle from the particle
%           pair forces.
%
% Notes:
%
% q, m, and alpha : arrays Nx1, or a function handles taking in one
%   parameter (the time), or a griddedInterpolant's taking in one parameter
%   (the time).
%
% dV : function handle or griddedInterpolant taking in
%   particle position and outputting the potential gradient at the
%   positions.
%
% dVint : function handle taking in the distance between two particles
%   and outputting the force.

% James Kapaldo

% Get all the inputs needed ==============================================
D = extraInputs.D; % dimension
dist = extraInputs.dist;
dist_arg = extraInputs.dist_arg;

% Particle properties
q = extraInputs.q; % particle charges (This charge should also include the square root of the coupling constant -- k or 1/(4*pi*eps0).)
m = extraInputs.m; % particle masses
alpha = extraInputs.alpha; % particle damping coefficients

if isa(q,'griddedInterpolant') || isa(q,'function_handle')
    q = q(t);
end
if isa(m,'griddedInterpolant') || isa(m,'function_handle')
    m = m(t);
end
if isa(alpha,'griddedInterpolant') || isa(alpha,'function_handle')
    alpha = alpha(t);
end

% Field properties and particle interactions.
dV = extraInputs.dV; % function handle taking in a location and giving dV/dx
dVint_fun = extraInputs.dVint; % function handle for particle interaction - takes in distance between two particles an outputs a scalar

% Indices for particle pairs.
pdistInds = extraInputs.pdistInds; % N*(N-1)/2 x 2 - indices linking each distance between two particles returned by pdist to the two particles.
N_to_NtNm1o2 = extraInputs.N_to_NtNm1o2;

% Number of inputs
[N,M] = size(y); % N = 2 * D * (number of particles), M = number of different input vectors

% Offset indices
offset = (0:2*D:N-1);

N = N/(2*D); % Number of particles

% Indices of position and velocity for each particle
rInds = (1:D)' + offset;
pInds = rInds + D;

r = y(rInds(:),:);
p = y(pInds(:),:);

% r is now a matrix that looks like
%  [ rx1_1, rx1_2, ..., rx1_M ;
%  [ ry1_1, ry1_2, ..., ry1_M ;
%  [   .  ,   .  , ...,   .
%  [ rD1_1, rD1_2, ..., rD1_M ;
%  [ rx2_1, rx2_2, ..., rx2_M ;
%  [ ry2_1, ry2_2, ..., ry2_M ;
%  [   .    .             .
%  [   .       .          .
%  [   .          .       .
%  [ rxN_1,      ...    rxN_M ;
%  [ ryN_1,      ...    ryN_M ;
%  [   .  ,      ...      .   ;
%  [ rDN_1, rDN_2, ..., rDN_M ];
%
% and p looks the same

% In order to calculate the potential at each nuclei position, we need the
% position vectors to be in the form
%
% [ x1,y1, ..., D1;            --> position vector of particle 1
%   x2,y2, ..., D2;            --> position vector of particle 2
%   x3,y3, ..., D3;
%    ... ;
%   xN*M,yN*M, ..., DN*M]      --> position vector of particle N*M
%
% (This is/(should be) the form expected by either a griddedInterpolant
% function or any other function handle created to give the potential
% forces.)

r = reshape(r,[D,N*M])';
dp = -dV(r);
% dp = zeros(N*M, D);
% for i = 1:D
%     dp(:,i) = -dV{i}(r);
% end

% In order to calculate the interaction between particles, we will reshape
% each set of particles to a page. Thus it will have N rows with two
% columns (x, y) and M pages
%
% Page 1
% [ rx1_1 ry1_1 ...;
% [ rx2_1 ry2_1 ...;
% [   .
% [   .
% [ rxN_1 ryN_1 ...];
%
%  .
%  .
%  .
%
% Page M
% [ rx1_M ry1_M ...;
% [ rx2_M ry2_M ...;
% [   .
% [   .
% [ rxN_M ryN_M ...];

r = permute( reshape( r', [D, N, M]), [2,1,3]); % NxDxM

% Convert dp into the same form.
dp = permute( reshape( dp', [D, N, M]), [2,1,3]); % NxDxM

% Also, put p into the same form
p = permute( reshape( p, [D, N, M] ), [2,1,3]); % NxDxM

% Get the distance between particles and save them as an N*(N-1)/2x1xM
% matrix

if M == 1
    % d = pdist(r);
    d = DN_pdistmex(r',dist,dist_arg)';
else
    d = zeros(N*(N-1)/2,1,M);

    for i = 1:M
        % d(:,1,i) = pdist(r(:,:,i));
        d(:,1,i) = DN_pdistmex(r(:,:,i)',dist,dist_arg);
    end
end
% Get the interaction force between the particles (which only depends on
% the distance between them)
dVint = dVint_fun(d); % N*(N-1)/2 x 1 x M

% Get the unit vectors between each pair

if M == 1
    r_hat = N_to_NtNm1o2*r; % N*(N-1)/2 x D
else
    r_hat = r(pdistInds(:,1),:,:) - r(pdistInds(:,2),:,:); % N*(N-1)/2 x D x M
end
r_hat = r_hat ./ sqrt(sum(r_hat.^2,2));

% Compute the force for each particle pair
if numel(q)==1
    forceOut = (r_hat .* dVint) .* q^2; % N*(N-1)/2 x D x M
else
    forceOut = (r_hat .* dVint) .* prod(q(pdistInds),2); % N*(N-1)/2 x D x M
end

% Accumulate the forces on each particle
if M == 1
    dp = dp + (-N_to_NtNm1o2')*forceOut;
else
    dpInt = zeros(N,D,M);
    for i = 1:M
        dpInt(:,:,i) = (-N_to_NtNm1o2')*forceOut(:,:,i);
    end
    dp = dp + dpInt;
end


% Add in the damping terms.
dp = dp - (alpha./m) .* p;

% Create the dr array.
dr = p ./ m;

% Reshape dr and dp
dr = reshape(permute(dr,[2,1,3]),D*N,M);
dp = reshape(permute(dp,[2,1,3]),D*N,M);

% Put dr and dp into dy, and then done!
dy = zeros(size(y));

dy(rInds(:),:) = dr;
dy(pInds(:),:) = dp;

end
