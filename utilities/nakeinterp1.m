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


% If no compiler is available to compile nakeinterp1.c, then we would need
% to include checks everywhere in the code to determine which interpolation
% function to use. Instead of doing that, here, we simply overload the
% nakeinterp1 function. If the .mex function is available, then it will be
% called; however, if it is not available, then the below code will run.
function yi = nakeinterp1(x, y, xi)

% 1D linear interpolation with linear extrapolation.
Y = griddedInterpolant(x,y);
yi = Y(xi);

end