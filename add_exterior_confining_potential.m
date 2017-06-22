function V = add_exterior_confining_potential(BW,V_interior,options)
% ADD_EXTERIOR_CONFINING_POTENTIAL
% Take the binary mask of a single object and the interior confining
% potential, pad them, and add in the exterior confining potential.
%
% V = add_exterior_confining_potential(BW,V_interior,options)
%
% Input parameters:
% BW         - The binary mask of single object
% V_interior - The confining potential inside the objects
% options    - Element of class seedPointOptions
%
% Output parameters:
% V - full confining potential
%
% Notes: We must pad the objects' binary mask image to have enough pixels
% around the edge to repel the particles back towards the object interior.
% This is necessary because the ODE solver does not take infinitesimally
% small time steps; thus, the PAD_SIZE should be larger than the distance
% that a particle could possibly jump in one time step.
%
% See also COMPUTESEEDPOINTS

% James Kapaldo

PAD_SIZE = options.Potential_Pad_Size;

% Process the mask -------------------------------------------------------

% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*[1 1]);
V = padarray(V_interior,PAD_SIZE*[1 1]);

% Set potential outside well to be strongly increasing with distance. Don't
% set to infinity because, when modeling, if a particle happened to go into
% the Inf region during a time step, then the solution will break.
V_out = (bwdist(BW_pad)+1);
V_out = V_out.^(V_out+1);
V(~BW_pad) = V_out(~BW_pad);

end
