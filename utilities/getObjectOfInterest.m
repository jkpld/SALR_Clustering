function varargout = getObjectOfInterest(options,varargin)
% GETOBJECTOFINTEREST  Return the elements in the cell arrays input in
% varargin that coorespond to options.Object_Of_Interest
%
% varargout = getObjectOfInterest(options,varargin)
% [A,B,C,...] = getObjectOfInterest(options,A_in,B_in,C_in,...)
%
% where A = A_in(options.Object_Of_Interest), ...

% James Kapaldo

if nargout > (nargin - 1)
    error('getObjectOfInterest:badInput','More outputs requested than input pixel lists given.')
end

if ~isempty(options.Object_Of_Interest)
    varargout = cell(1,numel(varargin));
    for i = 1:numel(varargin)
        varargout{i} = varargin{i}(options.Object_Of_Interest);
    end
else
    varargout = varargin;
end

end
