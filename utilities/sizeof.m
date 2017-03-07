function out = sizeof(in)

numclass = {'double'; 'single'; 'int8'; 'int16'; 'int32'; 'int64'; 'uint8'; 'uint16'; 'uint32'; 'uint64'};
numbytes = [8;4;1;2;4;8;1;2;4;8];

classIdx = strcmp(in,numclass);
if isempty(classIdx)
    error('sizeof:unknownClass','Input class, %s, is not a recognized class name', in)
end

out = numbytes(classIdx);

end