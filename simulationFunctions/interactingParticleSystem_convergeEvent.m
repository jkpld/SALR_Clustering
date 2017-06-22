function [value, isterminal, direction] = interactingParticleSystem_convergeEvent(~,y)
% INTERACTINGPARTICLESYSTEM_CONVERGEEVENT  Determine if the particle system
% has converged by checking if the mean speed is below some threshold.

% James Kapaldo

% Number of inputs (4 times the number of particles)
N = numel(y);

% Offset indices
offset = (0:4:N-1)';

% Indices of position and velocity for each particle
% rInds = [1,2] + offset;
pInds = [3,4] + offset;

% Positions and velocities of each particle
% r = y(rInds);
p = y(pInds);

value = mean(sqrt(sum(p.^2,2))) - 0.5e-2; % mean(speed) < THRESHOLD_SPEED
isterminal = 1;
direction = -1;

end