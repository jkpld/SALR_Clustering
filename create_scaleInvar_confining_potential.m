function [V, scaleFactor, overlapFactor] = create_scaleInvar_confining_potential(BW,options)
% CREATE_SCALEINVAR_CONFINING_POTENTIAL
% Take the binary mask of the image and create the base confining
% potential
%
% V = create_scaleInvar_confining_potential(BW,options)
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

PAD_SIZE = options.Potential_Padding_Size;
MAX_DT = options.Max_Distance_Transform;

% Process the mask -------------------------------------------------------
% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*ones(1,ndims(BW)));
BW_pad = logical(BW_pad);

% Compute bsae confining potential ---------------------------------------
% Compute the binary distance transform.
V = double(bwdist(~BW_pad));

% Extract object scale
object_scale = max(V(BW_pad));

% Extract the fraction of overlap in the object
% overlapFactor = get_overlap_factor(BW);
overlapFactor = 1;

% Scale confining potential if necessary, create external confining
% potential, and determing the scaleFactor.
if isnan(MAX_DT)
    scaleFactor = 1; % scale normalization factor
    overlapFactor = 1;
else
    % Normalize scale by scaling the distance transform
    V = 1 + (MAX_DT - 1) *(V - 1) / (object_scale - 1); % The +- 1's make sure that the distance transform goes from 1 to POTENTIAL_SCALE.
    scaleFactor = MAX_DT / object_scale; % scale normalization factor
end

% Add exterior confining potential
V_out = 1./(bwdist(BW_pad).^2 + 1);

% Insert the external confining potential
V(~BW_pad) = V_out(~BW_pad);

% Smooth the potential will small gaussian
if ismatrix(BW)
    V = imfilter(V,fspecial('gaussian',7,1),'replicate');
else
    % In higher dimensions we will filter in the frequency domain instead
    % of the spatial domain because it is faster.
    V = frequencyGaussianFilter(V,1,5,'replicate');
    
    % Sometimes this smoothing function can return very small magnitude
    % negative numbers where 0's should be. Since our confining potential
    % should be everywhere larger than zero, we can just set all negative
    % numbers to zero.
    V(V<0) = 0; 
end

% Re scale the potential to correct for the smoothing
if ~isnan(MAX_DT)
    V = MAX_DT * V / max(V(BW_pad));
end

% Invert to get final confining potential well.
V = 1./V;

end
