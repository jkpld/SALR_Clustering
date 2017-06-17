function varargout = getObjectOfInterest(options,varargin)
% GETOBJECTOFINTEREST  Return the options.Object_Of_Interest elements of
% the input arrays.
%
% varargout = getObjectOfInterest(options, varargin)
% [A, B, C] = getObjectOfInterest(options, A_in, B_in, C_in)
% A = A_in(options.Object_Of_Interest)
%
% Input parameters:
% options : An instance of class seedPointOptions
% varargin : Any number of arrays or cell arrays
%
% Output parameters:
% varargout : The options.Object_Of_Interest element of each input.

% James Kapaldo

if nargout > (nargin - 1)
    error('getObjectOfInterest:badInput','More outputs requested than input pixel lists given.')
end

if ~isempty(options.Object_Of_Interest)
    varargout = cell(1,numel(varargin));
    for i = 1:numel(varargin)
        if isempty(varargin{i})
            varargout{i} = [];
        else
            varargout{i} = varargin{i}(options.Object_Of_Interest);
        end
    end
else
    varargout = varargin;
end

end
