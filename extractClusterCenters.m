function clusterCenters = extractClusterCenters(r_final,options)
% EXTRACTCLUSTERCENTERS Extract the center location of each cluster formed
% by the particles.
%
% Input parameters:
% r_final - The final locations of the particles, as returned by
%           modelParticleDynamics
% options - Element of class seedPointOptions
%
% Output parameters:
% clusterCenters - N x size(r_final,2) array where each row is a cluster
%                  center.
%
% See also MODELPARTICLEDYNAMICS

% James Kapaldo
% 2017-01-20

% Number of particles
N = size(r_final,1);
D = size(r_final,2);

% Distance between each pair of particles 
d = pdist(r_final);

% Particle pairs that are connected to each other
cluster = (d * options.Scale_Factor) < (0.7*options.Potential_Minimum_Location + 0.3*options.Potential_Extent);
pdistInds = getPdistInds(N);
linIdx = pdistInds(cluster,1) + (pdistInds(cluster,2)-1)*N;

% Symmetric particle adjacency matrix
A = eye(N,'logical');
A(linIdx) = true;
A = A | A';

% Permute to block diagonal form and extract the size of each block
[p,~,r] = dmperm(A);

% Compute the center of each cluster.
clusterCenters = zeros(numel(r)-1,D);
for i = 1:numel(r)-1
    idx = r(i):r(i+1)-1;
    clusterCenters(i,:) = mean( r_final( p(idx) ,:), 1);
end

end
