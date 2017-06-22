function [hx,hy] = cart2hex(cp,a,b)
% CART2HEX  Convert from cartesian points to hexagon center points.
% Hexagons will have a lattice constant of A and a aspect ratio of B.
% Output will be an array of the same size as cp with the hexagonal centers
% of each point. Hexagon center 0,0 is at cartesian point 0,0.
%
% [hx,hy] = cart2hex(cp,a,b)
%
% See also HEX2CART, HEXPLOT

% Ref: 
% http://playtechs.blogspot.co.uk/2007/04/hex-grids.html
% https://gist.github.com/egradman/583180


R = a/2;
S = b*2*R/sqrt(3);
% T = S/2;

% Tx = [0 R; -S, S/2];
iTx = [1/(2*R), -1/S; 1/R, 0];

% Ty = [R, -R; S/2, S/2];
iTy = [1/(2*R), 1/S; -1/(2*R), 1/S];


hx = [floor((sum(floor(iTx*cp'),1)' + 2)/3), floor((sum(floor(iTy*cp'),1)' + 2)/3)];

if nargout > 1
    hy = hx(:,2);
    hx = hx(:,1);
end

