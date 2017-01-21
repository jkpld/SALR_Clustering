function clusterCenters = extractClusterCenters(r_final,objSize,options)
% EXTRACTCLUSTERCENTERS Extract the center location of each cluster formed
% by the particles.
%
% Input parameters:
% r_final - The final locations of the particles, as returned by
%           modelParticleDynamics
% objSize - The size of the object mask image
% options - Element of class seedPointOptions
%
% Output parameters:
% clusterCenters - N x size(r_final,2) array where each row is a cluster
%                  center.
%
% See also MODELPARTICLEDYNAMICS

% James Kapaldo
% 2017-01-20

% Linear indices of particle locations
rInd = r_final(:,1) + (r_final(:,2)-1) * objSize(1);

rInd(~isfinite(rInd)) = [];

% Search for the clusters by creating an image with the particles as 1's,
% then dilate and find the centroids of the connected components.

mask = false(objSize);
mask(rInd) = 1;
mask = imdilate(mask,strel('disk',ceil(options.Potential_Minimum_Location)));

CC = bwconncomp(mask);
props = regionprops(CC,'Centroid'); % Perhaps should replace this with faster code for getting centroids.

% Combine all of the centroids into one array and flip the xy components as
% regionprops changes the order.
clusterCenters = fliplr(cat(1,props.Centroid));

end
