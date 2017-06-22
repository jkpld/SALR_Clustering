function [cpx,cpy] = hex2cart(hx,a,b)
% HEX2CART  Convert from hexagon center points to cartesian points.
% The hexagons will have a lattice constant of A and a aspect ratio of B.
% Output will be an array of the same size as hp with the cartesian
% coordiantes of each hexagon center. Hexagon center 0,0 is at cartesian
% point 0,0.
%
% [cpx,cpy] = hex2cart(hx,a,b)
%
% See also CART2HEX, HEXPLOT

% Ref: 
% http://playtechs.blogspot.co.uk/2007/04/hex-grids.html
% https://gist.github.com/egradman/583180


R = a/2;
S = b*2*R/sqrt(3);
T = S/2;

% Tx = [0 R; -S, S/2];
% iTx = [1/(2*R), -1/S; 1/R, 0];

% Ty = [R, -R; S/2, S/2];
% iTy = [1/(2*R), 1/S; -1/(2*R), 1/S];


% hx = [floor((sum(floor(iTx*cp'),1)' + 2)/3), floor((sum(floor(iTy*cp'),1)' + 2)/3)];

cpx = [hx(:,1)*(2*R)+hx(:,2)*R, hx(:,2)*(S+T)];

if nargout > 1
    cpy = cpx(:,2);
    cpx = cpx(:,1);
end