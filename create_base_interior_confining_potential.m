function [V, scaleFactor] = create_base_interior_confining_potential(BW,options)
% CREATE_BASE_INTERIOR_CONFINING_POTENTIAL
% Take the binary mask of the image and create the base inteior confining
% potential
%
% V = create_base_interior_confining_potential(BW,options)
%
% Input parameters:
% BW        - The binary mask
% options   - Element of class seedPointOptions
%
% Output parameters:
% V - Base interior confining potential
%
% See also COMPUTESEEDPOINTS

% James Kapaldo

PAD_SIZE = 1;
POTENTIAL_SCALE = options.Potential_Scale;

% Process the mask -------------------------------------------------------

% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*[1 1]);

% Compute bsae confining potential ---------------------------------------

% Get potential inside well. Compute the binary distance transform, smooth
% it out, and invert it to get 1/r.
V = double(bwdist(~BW_pad));
% V = imfilter(V,fspecial('gaussian',7,1)); % Move this filter to the add_external_confining_potential function

% Extract object scale 
object_scale = max(V(:));

% Enusre that the first pixel in from the boundary has a value of 1 at the
% end.
V = V - 1; % min(V(BW_pad)); % if the above gaussian filtering is used, then use the min(V(BW_pad))) instead of 1.
V(V<0) = 0;

if isnan(POTENTIAL_SCALE)
    V = V + 1;
    scaleFactor = 1; % scale normalization factor
else
    % Normalize scale by scaling the distance transform
    V = (POTENTIAL_SCALE-1)*V./(object_scale - 1) + 1;    
    scaleFactor = POTENTIAL_SCALE / object_scale; % scale normalization factor
end

V = 1./V;

% Remove padding
V = V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);

end
