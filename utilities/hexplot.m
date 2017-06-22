function ph = hexplot(cartCents,a,b,varargin)
% HEXPLOT  Create a hexagonal scatter plot. 
% Hexagons with lattice constant A and aspect ratio B will be plotted at
% each cartesian center given. Hexagon coloring and relative size can be
% controled through parameters in options.
%
% ph = hexplot(cartCents,a,b,varargin)
%
% Input parameters:
% 
% cartCents - Nx2 array where each row gives the x,y location of a hexagon
%             center. (This data would normally be returned by hex2cart().)
%
% a - lattice constant of hexagons.
%
% b - aspect ratio of hexagons. (This will streach out hexagons along the
%     second dimenstion (y).)
%
% Optional parameters:
%
% colorData - Nx1 array giving the color of each hexagon in terms of the
%             colormap. (This will be used with the FaceVertexCData
%             property of the patches.) Instead of giving an Nx1 array you
%             can give the string 'sameAsSize' to use the same data given
%             for the size. Use 'none' to to have no color data (in this
%             case all hexagons will be black). Default is 'none'.
%
% sizeData - Nx1 array with values in the range of (0 1] giving the
%            relative sizing of each hexagon. Instead of giving an Nx1
%            array you can give the string 'sameAsColor' to use the same
%            data given for the size. Use 'none' to have no size data.
%            Default is 'none'.
%
% colorScale - Scale applied to the color data. Options are 'linear' or
%              'log'. Default is 'linear'.
%
% sizeScale - Scale applied to the size data. Options are 'linear', 'log', 
%             or 'none'. Default is 'none' if sizeData is not given and
%             'linear' if sizeData is given. Note, if 'none' is specified
%             then sizeData will be ignored.
%
% maxHexSize - A scalar from (0,1]. If maxHexSize is 1, then full size
%              hexagons (either hexagons with sizeData equal to 1 or when
%              no sizeData is given) will touch each other. If maxHexSize
%              is less than 1, then there will always be space between each
%              of the hexagons. Default is 0.9;
%
% orientation - The orientation of the hexagons in degrees. Default is 0,
%               which has the flat edges parallel along y.
%
% parent - Parent graphics object to put the patch object in. Default is
%          to create a new figure.
%
% Output parameters:
%
% p - handle to the patch object created.
%
% See also CART2HEX, HEX2CART

% JamesKapaldo
% 2016-10-10

% Get input values
ip = inputParser;
ip.FunctionName = 'hexplot';
addParameter(ip,'colorData','none',@(t) (size(t,1) == size(cartCents,1)) || any(strcmp(t,{'none','sameAsSize'})));
addParameter(ip,'sizeData','none',@(t) ((size(t,1) == size(cartCents,1)) && (size(t,2) == 1)) || any(strcmp(t,{'none','sameAsColor'})));
addParameter(ip,'colorScale','linear',@(t) any(strcmp(t,{'linear','log'})))
addParameter(ip,'sizeScale','default',@(t) any(strcmp(t,{'linear','log','none','default'})))
addParameter(ip,'maxHexSize',0.9,@(t) t > 0 && t<=1);
addParameter(ip,'orientation',0,@(t) isfinite(t) && numel(t)==1);
addParameter(ip,'parent',1)


parse(ip,varargin{:})

colorData = ip.Results.colorData;
sizeData = ip.Results.sizeData;
colorScale = ip.Results.colorScale;
sizeScale = ip.Results.sizeScale;
maxHexSize = ip.Results.maxHexSize;
orientation = ip.Results.orientation;
parent = ip.Results.parent;

if strcmp(colorData,'sameAsSize') && strcmp(sizeData,'sameAsColor')
    error('hexplot:badSizeColorData','You cannot use ''sameAsSize'' and ''sameAsColor''. You must give data for at least one of them or use ''none''.');
end

if strcmp(sizeData,'none')
    sizeData = [];
end

if strcmp(colorData,'none')
    colorData = [];
end

if strcmp(colorData,'sameAsSize')
    colorData = sizeData;
end

if strcmp(sizeData,'sameAsColor')
    sizeData = colorData;
end

if any(strcmp(sizeScale,{'linear','log'})) && isempty(sizeData)
    warning('hexplot:sizeScaleSetWithNoSizeData','sizeData not given, sizeScale being ignored.')
    sizeScale = 'none';
end

if ~isempty(sizeData) && strcmp(sizeScale,'default')
    sizeScale = 'linear';
end

R = a/2;
S = b*2*R/sqrt(3);

hexVerts = [ 0,  S;
            -R,  S/2;
            -R, -S/2;
             0, -S;
             R, -S/2;
             R,  S/2];

R = [cosd(orientation), -sind(orientation); sind(orientation), cosd(orientation)];
hexVerts = (R*hexVerts')';
         
if ~isempty(sizeData) && ~strcmp(sizeScale,'none')
    switch sizeScale
        case 'linear'
            nS = sizeData;
        case 'log'
            nS = max(log(sizeData),0);
    end
    nS = nS/max(nS);
    hexVerts = hexVerts .* permute(nS,[3,2,1]);
end     

N = size(cartCents,1);

hexVerts = maxHexSize * hexVerts + permute(cartCents,[3,2,1]);
hexVerts = reshape(permute(hexVerts,[2,1,3]),2,6*N)';


hexFaces = reshape(1:6*N,6,N)';

if parent == 1
    fig = figure;
    parent = axes('parent',fig);
end


if isempty(colorData)
    p = patch('Vertices',hexVerts,'Faces',hexFaces,'FaceColor','k','EdgeColor','none','Parent',parent);
else
    p = patch('Vertices',hexVerts,'Faces',hexFaces,'FaceColor','flat','EdgeColor','flat','Parent',parent);
    switch colorScale
        case 'linear'
            n = colorData;
        case 'log'
            n = max(log(colorData),0);
    end
    p.FaceVertexCData = kron(n,ones(6,1));
end

if nargout > 0
    ph = p;
end

end