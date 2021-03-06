function out = sizeof(in)
% SIZEOF Get the size in bytes of the input numeric class.
%
% out = sizeof(in)
%
% Input parameters:
% in : The name of a numeric type
%
% Output parameters:
% out : The size of the input type in bytes.

% James Kapaldo

numclass = {'double'; 'single'; 'int8'; 'int16'; 'int32'; 'int64'; 'uint8'; 'uint16'; 'uint32'; 'uint64'};
numbytes = [8;4;1;2;4;8;1;2;4;8];

classIdx = strcmp(in,numclass);
if isempty(classIdx)
    error('sizeof:unknownClass','Input class, %s, is not a recognized class name', in)
end

out = numbytes(classIdx);

end


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
