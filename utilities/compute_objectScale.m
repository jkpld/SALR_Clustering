function objectScale = compute_objectScale(BW, pixelList)
% COMPUTE_OBJECTSCALE Compute the maximum distance transform value for
% each object in a binary mask.
%
% objectScale = compute_objectScale(BW, pixelList)

% James Kapaldo

PAD_SIZE = 1;

% Pad the binary mask with with 0's so the distance transform is correct
% for objects on the edge.

BW_pad = padarray(BW,PAD_SIZE*[1 1]);
BW_pad = logical(BW_pad);

% Compute distance transform
DT = double(bwdist(~BW_pad));

% Remove padding
DT = DT(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE);

% Get the object scale, which is just the maximum distance transform value
% for each object.
objectScale = cellfun(@(x) max(DT(x)), pixelList);
objectScale = objectScale(:);

% L = bwlabel(BW);
% toRemove = L==0;
% DT(toRemove) = [];
% L(toRemove) = [];
% objectScale = accumarray(L,DT,[max(L),1],@max);


end
