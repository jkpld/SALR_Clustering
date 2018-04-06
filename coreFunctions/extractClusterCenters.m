function [clusterCenters, clusterSize] = extractClusterCenters(r_final,options)
% EXTRACTCLUSTERCENTERS Extract the center and the number of particles in
% particle cluster.
%
% [clusterCenters, clusterSize] = extractClusterCenters(r_final,options)
%
% Input parameters:
% r_final : The final locations of the particles, as returned by
%           modelParticleDynamics
% options : Element of class seedPointOptions
%
% Output parameters:
% clusterCenters : N x D array where each row is a cluster center. (D is 
%                  the dimension).
% clusterSize : N x 1 array where the i'th element gives the number of
%               particles in the i'th cluster
%
% See also MODELPARTICLEDYNAMICS

% James Kapaldo
% 2017-01-20

% Parameters
ra = options.ScaleInvarient_Potential_Extent;
r0 = options.ScaleInvarient_Potential_Minimum_Location;
if isnan(ra)
    ra = options.Potential_Parameters(3);
    r0 = options.Potential_Parameters(2);
end


dist = options.dist;
dist_arg = options.dist_arg;

% Number of particles
N = size(r_final,1);
D = size(r_final,2);

% Distance between each pair of particles (in data units)
d = pdist(r_final,dist,dist_arg);

% Particle pairs that are connected to each other
cluster = d  < (0.7*r0 + 0.3*ra);
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
clusterSize = zeros(numel(r)-1,1);
for i = 1:numel(r)-1
    idx = r(i):r(i+1)-1;
    clusterSize(i) = length(idx);
    clusterCenters(i,:) = mean( r_final( p(idx) ,:), 1);
end

end


%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
