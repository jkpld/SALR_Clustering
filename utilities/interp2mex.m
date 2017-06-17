% INTERP2MEX Fast 2-D bilinear interpolation
%
% Zi = interp2mex(Z, Xi, Yi)
%
% Interpolates to find Zi, the values for the underlying 2-D function Z at
% the points in matrices Xi and Yi. Out of range values are interpolated to
% nearest neighbor. Only input matrices of class double are supported.
% Other classes will give unpredictable results.
%
% This function assumes that the spacing between all points is 1.
%
% The code is built starting from the function posted on the fileexchange
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/68708
% by Jøger Hansegård

% James Kapaldo
