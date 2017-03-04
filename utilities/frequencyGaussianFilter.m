function A = frequencyGaussianFilter(A, sigma, hsize, padding)


% This function is a modified form of the function by the same name in
% imgaussfilt3

dtype = class(A);
outSize = size(A);

pdsz = floor(hsize/2) * ones(1, ndims(A));
A = padarray(A,pdsz,padding,'both');

h = ndGaussianFilter(ndims(A),sigma,hsize);

% cast to double to preserve precision unless single
if ~isfloat(A)
    A = double(A);
end

fftSize = size(A);
A = ifftn( fftn(A, fftSize) .* fftn(h, fftSize), 'symmetric' );

if ~strcmp(dtype,class(A))
    A = cast(A, dtype);
end

start = 1 + size(A)-outSize;
stop = start + outSize - 1;

inds = cell(1,ndims(A));
for i = 1:ndims(A)
    inds{i} = start(i):stop(i);
end

A = A(inds{:});

end
