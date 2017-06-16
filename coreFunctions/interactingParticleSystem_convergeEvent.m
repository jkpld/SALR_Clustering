function [value, isterminal, direction] = interactingParticleSystem_convergeEvent(t,y,m,hstry,D)
% INTERACTINGPARTICLESYSTEM_CONVERGEEVENT  Determine if the particle system
% has converged by checking if the mean particle speed is below some
% threshold.
%
% [value, isterminal, direction] = interactingParticleSystem_convergeEvent(t,y,m,hstry,D)

% James Kapaldo

% Number of inputs (2D times the number of particles)
N = numel(y);

% Offset indices
offset = (0:2*D:N-1)';

N = N/(2*D); % Number of particles

% Indices of position and velocity for each particle
% rInds = [1,2] + offset;
pInds = (1:D) + D + offset;

% Positions and velocities of each particle
% r = y(rInds); % N x 2
p = y(pInds);

% ================  NaN error quit ==================================
if any(isnan(y))
    value(2) = -1;
else
    value(2) = 1;
end
isterminal(2) = 1;
direction(2) = -1;

% ================  Get history if full =============================
W = hstry.windowSize;

if W > 0

    if  hstry.validEntries >= W
        pInds = pInds';
        p_hist = hstry.solHistory(pInds(:),:); % 2N x windowSize
        p_hist = permute( reshape( p_hist, [D, N, W] ), [2,1,3]); % Nx2xM
    end

    % ================  AVERAGE POSITION STABLE  ========================
    % if  hstry.iterationNumber < W
    %     value(2) = 1;
    % else
    %     d = sum((r_hist - r).^2,2);
    %     value(2) = mean(d(:)) - 0.01;
    % end
    % isterminal(2) = 1;
    % direction(2) = -1;

    % ===============  AVERAGE VELOCITY STABLE =========================
    if  hstry.validEntries < W
        value(1) = 1;
    else
        v = sum(mean(cat(3,p_hist,p),3).^2,2); % average speed of each particle over history
        value(1) = mean(v(:))/m - 0.2e-3; % mean(speed) < THRESHOLD_SPEED
    end
    isterminal(1) = 1;
    direction(1) = -1;

    hstry.add(t,y);
else
    % ================  SPEED MAGNITUDE CONVERGE EVENT ==================
    value(1) = mean(sqrt(sum(p.^2,2)))./m - 0.5e-2; % mean(speed) < THRESHOLD_SPEED
    isterminal(1) = 1;
    direction(1) = -1;
end


end
