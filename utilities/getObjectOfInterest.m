function varargout = getObjectOfInterest(Object_Of_Interest,varargin)
% GETOBJECTOFINTEREST  Return the Object_Of_Interest element of the input
% arrays.
%
% varargout = getObjectOfInterest(Object_Of_Interest, varargin)
% [A, B, C] = getObjectOfInterest(Object_Of_Interest, A_in, B_in, C_in)
% A = A_in(Object_Of_Interest)
%
% Input parameters:
% Object_Of_Interest : Index of the object of interest
% varargin : Any number of arrays or cell arrays
%
% Output parameters:
% varargout : The options.Object_Of_Interest element of each input.

% James Kapaldo

if nargout > (nargin - 1)
    error('getObjectOfInterest:badInput','More outputs requested than input pixel lists given.')
end

if ~isempty(Object_Of_Interest)
    varargout = cell(1,numel(varargin));
    for i = 1:numel(varargin)
        if isempty(varargin{i})
            varargout{i} = [];
        else
            varargout{i} = varargin{i}(Object_Of_Interest);
        end
    end
else
    varargout = varargin;
end

end
