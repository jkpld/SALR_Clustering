function [V, scaleFactor,object_scale] = create_scaleInvar_confining_potential(BW,options)
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
POTENTIAL_SCALE = options.Potential_Scale;

% Process the mask -------------------------------------------------------
% Pad mask with zeros so that we have a strong potential all around the
% object.
BW_pad = padarray(BW,PAD_SIZE*[1 1]);
BW_pad = logical(BW_pad);

% Compute bsae confining potential ---------------------------------------
% Compute the binary distance transform.
V = double(bwdist(~BW_pad));

% Extract object scale
object_scale = max(V(BW_pad));
% object_scale = prctile(V(logical(BW_pad)),90);
% object_scale = mean(V(  (V==imdilate(V,strel('disk',2))) & (V>0) ));

% Scale confining potential if necessary, create external confining
% potential, and determing the scaleFactor.
if isnan(POTENTIAL_SCALE)
    scaleFactor = 1; % scale normalization factor

    % Add exterior confining potential
    V_out = 1./(bwdist(BW_pad).^2 + 1);
else
    % Normalize scale by scaling the distance transform
    V = POTENTIAL_SCALE*(V/object_scale);
    scaleFactor = POTENTIAL_SCALE / object_scale; % scale normalization factor

    % Add exterior confining potential
    V_out = 1./(bwdist(BW_pad).^2 + object_scale/POTENTIAL_SCALE);
end

% Insert the external confining potential
V(~BW_pad) = V_out(~BW_pad);

% Smooth the potential will small gaussian
V = imfilter(V,fspecial('gaussian',7,1));

% Re scale the potential to correct for the smoothing
if ~isnan(POTENTIAL_SCALE)
    V = POTENTIAL_SCALE * V / max(V(BW_pad));
end

% Invert to get final confining potential well.
V = 1./V;

end
