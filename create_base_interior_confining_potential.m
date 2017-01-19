function V = create_base_interior_confining_potential(BW)
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

% Process the mask -------------------------------------------------------

% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*[1 1]);

% Compute bsae confining potential ---------------------------------------

% Get potential inside well. Compute the binary distance transform, smooth
% it out, and invert it to get 1/r.
V = double(bwdist(~BW_pad));
V = imfilter(V,fspecial('gaussian',7,1));
V = 1./V;

% Remove padding
V = V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);

end
