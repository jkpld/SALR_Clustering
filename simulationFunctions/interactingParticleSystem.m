function dy = interactingParticleSystem(t,y,extraInputs)
% INTERACTINGPARTICLESYSTEM  Given particle position (r) and momentum (p)
% compute dr/dt and dp/dt.
%
% The parameter vector y should have the form
%
% y  = [y1 | y2 | y3 | ... | yn];
% y1 = [ rx1, ry1, vx1, vy1, rx2, ry2, vx2, vy2, ... , rxN, ryN, vxN, vyN ]'
% 
% Where N is the number of particles and where n > 1 if the ode solver is
% trying to compute multiple points at once.
%
% ExtraInputs - structure with the fields
%       q - particle charges
%       m - particle masses
%       alpha - particle damping coefficients
%       dVx - x derivative of potential
%       dVy - y derivative of potential
%       dVind - derivative of particle interaction potential
%       pdistInds - indices of particle pairs
%       N_to_NtNm1o2 - ( read as : N to N*(N-1)/2 ) sparse matrix that can
%       create particle pair vectors from particle positions and accumulate
%       all of the forces on a particle from the particle pair forces.
%
% Notes: 
% 
% q, m, and alpha  :  arrays Nx1, a function handel taking in one paramter,
% the time, or a griddedInterpolant taking in one paramter, the time.
%
% dVx, dVy  :  function handle or griddedInterpolant taking in particle
% position and outputing dVx or dVy
%
% dVint  :  function handel taking in the distance between two particles
% and outputing the force.

% James Kapaldo

% Get all the inputs needed ==============================================
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
dVx = extraInputs.dVx; % function handle taking in a location and giving dV/dx
dVy = extraInputs.dVy; % function handle taking in a location and giving
dVint_fun = extraInputs.dVint; % function handle for particle intercation - takes in distance between two particles an outputs a scalar

% Indices for particle pairs.
pdistInds = extraInputs.pdistInds; % N*(N-1)/2 x 2 - indices linking each distance between two partices returned by pdist to the two particles.
N_to_NtNm1o2 = extraInputs.N_to_NtNm1o2;

% Number of inputs
[N,M] = size(y); % N = 4 * (number of particles), M = number of different input vectors

% Offset indices
offset = (0:4:N-1);

N = N/4; % Number of particles

% Indices of position and velocity for each particle
rInds = [1;2] + offset;
pInds = [3;4] + offset;

r = y(rInds(:),:);
p = y(pInds(:),:);

% r is now a matrix that looks like
%  [ rx1_1, rx1_2, ..., rx1_M ;
%  [ ry1_1, ry1_2, ..., ry1_M ;
%  [ rx2_1, rx2_2, ..., rx2_M ;
%  [ ry2_1, ry2_2, ..., ry2_M ;
%  [   .    .             .
%  [   .       .          .
%  [   .          .       .
%  [ rxN_1,      ...    rxN_M ]
%  [ ryN_1,      ...    ryN_M ]
%
% and p looks the same

% In order to calculate the potential at each nuclei position, we need the
% position vectors to be in the form 
%
% [ x1,y1; 
%   x2,y2; 
%   x3,y3; 
%    ... ; 
%   xN*M,yN*M]
%
% (This is/(should be) the form expected by either a griddedInterpolant
% function or any other function handle cretaed to give the potentials.)

r = reshape(r,[2,N*M])';
dp = -[dVx(r),dVy(r)];

% In order to calculate the interaction between particles, we will reshape
% each set of particles to a page. Thus it will have N rows with two
% columns (x, y) and M pages
% 
% Page 1
% [ rx1_1 ry1_1 ;
% [ rx2_1 ry2_1 ;
% [   .
% [   .
% [ rxN_1 ryN_1 ];
%
%  .
%  .
%  .
% 
% Page M
% [ rx1_M ry1_M ;
% [ rx2_M ry2_M ;
% [   .
% [   .
% [ rxN_M ryN_M ];

r = permute( reshape( r', [2, N, M]), [2,1,3]); % Nx2xM

% Convert dp into the same form.
dp = permute( reshape( dp', [2, N, M]), [2,1,3]); % Nx2xM

% Also, put p into the same form
p = permute( reshape( p, [2, N, M] ), [2,1,3]); % Nx2xM

% Get the distance between particles and save them as an N*(N-1)/2x1xM
% matrix

if M == 1
    % D = pdist(r);
    D = DN_pdistmex(r','euc',[])';
else
    D = zeros(N*(N-1)/2,1,M);

    for i = 1:M
        % D(:,1,i) = pdist(r(:,:,i));
        D(:,1,i) = DN_pdistmex(r(:,:,i)','euc',[]);
    end
end
% Get the interaction force between the particles (which only depends on
% the distance between them)
dVint = dVint_fun(D); % N*(N-1)/2 x 1 x M

% Get the unit vectors between each pair

if M == 1
    r_hat = N_to_NtNm1o2*r; % N*(N-1)/2 x 2
else
    r_hat = r(pdistInds(:,1),:,:) - r(pdistInds(:,2),:,:); % N*(N-1)/2 x 2 x M
end
r_hat = r_hat ./ sqrt(sum(r_hat.^2,2));

% Compute the force for each particle pair
if numel(q)==1
    forceOut = (r_hat .* dVint) .* q^2; % N*(N-1)/2 x 2 x M
else
    forceOut = (r_hat .* dVint) .* prod(q(pdistInds),2); % N*(N-1)/2 x 2 x M
end

% Accumulate the forces on each particle
if M == 1
    dp = dp + (-N_to_NtNm1o2')*forceOut;
else
    dpInt = zeros(N,2,M);
    for i = 1:M
        dpInt(:,:,i) = (-N_to_NtNm1o2')*forceOut(:,:,i);
    end
    dp = dp + dpInt;
end


% Add in the damping terms.
dp = dp - (alpha./m) .* p;

% Create the dr array.
dr = p ./ m;%  - (alpha./m) .* p;

% Limit the maximum particle speed for stability
dr(dr<-1) = -1;
dr(dr>1) = 1;

% Reshape dr and dp
dr = reshape(permute(dr,[2,1,3]),2*N,M);
dp = reshape(permute(dp,[2,1,3]),2*N,M);

% Put dr and dp into dy, and then done!
dy = zeros(size(y));

dy(rInds(:),:) = dr;
dy(pInds(:),:) = dp;

end





