function A = removePadding(A,pad_size)
% REMOVEPADDING Remove padding from and N-D matrix.
%
% A = removePadding(A, pad_size)
%
% Input parameters:
% A : input N-D matrix
% pad_size : the amount of padding to removed (symmetrically)
%
% Output parameters:
% A : the input matrix with the padding removed.

% James Kapaldo

sz = size(A);
N = length(sz);

if N ~= numel(pad_size)
    pad_size = pad_size(1)*ones(1,N);
end

inds = cell(1,N);

for i = 1:ndims(A)
    inds{i} = (pad_size(i)+1) : (sz(i)-pad_size(i));
end

A = A(inds{:});

end


%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
