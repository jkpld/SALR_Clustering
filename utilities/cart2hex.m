function [hx,hy] = cart2hex(cp,a,b)
% CART2HEX Convert from cartesian points to hexagon center points.
%
% [hx,hy] = cart2hex(cp,a,b)
%
% Input parameters:
% cp : Cartesian points, Nx2
% a : Hexagon lattice constant
% b : Hexagon aspect ratio
%
% Output parameters:
% hx, hy : Hexagon center points, Nx2
%
% Note:
% Hexagon center 0,0 is at Cartesian point 0,0.
%
% See also HEX2CART, HEXPLOT

% James Kapaldo

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


%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
