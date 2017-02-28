function varargout = createObjectImages(pixelList,N,varargin)
% CREATEOBJECTIMAGES  Create the mask, edge, and intensity image for the
% object with pixels given by pixelList from an image with N number
% of rows.

% James Kapaldo

if nargout > (nargin - 2)
    error('createObjectImages:badInput','More outputs requested than input pixel lists given.')
end

% Get linear indices of object
j = rem(pixelList-1, N) + 1;
k = (pixelList - j)/N + 1;

r = [j,k];

% remove top left corner
minR = min(r)-1;
r = r-minR;

% Get the boundary size
maxR = max(r,[],1,'omitnan');

% Get linear indices for small image
inds = r(:,1) + (r(:,2)-1)*maxR(1);

% Initialize output
varargout = cell(1,numel(varargin));

% Assign output
for i = 1:numel(varargin)
    if isempty(varargin{i})
        varargout{i} = [];
    else
        tmp = zeros(maxR);
        tmp(inds) = varargin{i};
        varargout{i} = tmp;
    end
end
end
% createObjectImages changes log
% 2016-10-17
% 2016-10-18/28 : undocumented changes
% 2016-10-29 : changed name to createObjectImages, renamed some inputs.
% 2017-01-19 : rewrote to include arbitrary number of pixel value lists
% using varargin
