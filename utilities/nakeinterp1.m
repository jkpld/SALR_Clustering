% NAKEINTERP1 Fast 1D interpolation using a dichotomy search of indices.
%
% idx = nakeinterp1(x, y, xi);
%
% where x, y and xi are double column vectors
% x must be sorted in ascending order; x and y have the same length
%
% NO ARGUMENT CHECKING
% Compile:
% mex -O -v nakeinterp1.c
% Author: Bruno Luong
% Original: 19/Feb/2009
