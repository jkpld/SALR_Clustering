function h = ndGaussianFilter(N,sigma,hsize)
% NDGAUSSIANFILTER Create an N-D Gaussian filter.
%
% h = ndGaussianFilter(N,sigma,hsize)
%
% Input parameters:
% N : number of dimensions
% sigma : Gaussian sigma
% hsize : filter size
%
% Output parameters:
% h : N-D array with hsize elements along each dimension

% James Kapaldo

% Compute filter size from sigma if filter size is not given
if nargin < 3
    hsize = 2*ceil(2*sigma) + 1;
end

% Filter radius
r = (hsize-1)/2;

% Create coordinate array and get length
d = (-r : r)';
n = length(d);

% initialize filter
h = zeros(n*ones(1,N));

% Create argument of exponential
d = d.^2 / (2* sigma^2);
order = 1:N;
for i = 1:N
    h = h + permute(d,circshift(order,-i+1));
end

h = exp(-h);

% Suppress near-zero components
h(h<eps*max(h(:))) = 0;

% Normalize
sumH = sum(h(:));
if sumH ~=0
    h = h./sumH;
end


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
