function h = ndGaussianFilter(N,sigma,hsize)

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