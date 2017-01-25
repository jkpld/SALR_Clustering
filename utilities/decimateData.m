function dcmt = decimateData(x,y,z,varargin)
% DECIMATEDATA  Bin 2D spatial data and return a single value for each bin.
% 
% dcmt = decimateData(x,y,z) will take spatial data -scatter data- (x,y,z)
% and decimate it to a single xy rectangular grid with grid spacing
% binSize. Specifically, the x-y data will be binned and the z data in each
% bin will be operated on to return the median value. The output is a
% structure array with fields, X, Y, and Z.
%
% Input parameters:
% x - Spatial x data
% y - Spatial y data
% z - data.
%
% dcmt = decimateData(x,y,z, ...) Custumize the binSize, the type of grid,
% the method of reduction, and choose if you want to produce a plot. These
% options can be given using parameter/value pairs. The parameters are
% described below.
%
% Optional parameters:
%
% binSize - Size of the spatial bin to bin x and y data. Default is 0.25.
%
% gridType - The type of the grid uesed. Available grids are 'rectangular'
% and 'hexagonal'. In 'rectangular', binSize(1) and binSize(2) give the
% side lengths of the rectangle. In 'hexagonal', binSize(1) gives the
% lattice constant of the hexagonal grid and binSize(2) gives the aspect
% ratio of the hexagons. Ex. if binSize(1) = 0.5 then hexagons with a
% lattice constant of 0.5 are created. After creation, the hexagons are all
% streatched along the second dimension (y) by a factor of binSize(2).
% Default is 'rectangular'.
%
% reductionMethod - The method of reducing all data in one spatial bin to a
% single point. Available methods are 'mean', 'median', 'mode',
% 'percentile', 'std', 'mad', 'density', or a custom function handle that
% takes in an array and outputs a single number. If mode is selected, then
% the data is first rounded to the nearest MAD/50, where MAD ('mad') is the
% median abosolute deviation. If percentile is selected then,
% prctile(z,PRCNT), will be used. You can set PRCNT with the optional
% parameter 'percentValue'. If density is selected then the z array is not
% used and the density of the x, y data is computed. Default is median.
%
% cleanHexagonData - boolian, only applicable if 'gridType' is 'hexagonal'.
%                    If true, then hexagons without any data points inside
%                    will be removed. Default is false.
%                     
%                    This option is only given for hexagonal, since the
%                    data will need to be interpolated with a scattered
%                    interpolant anyway. In rectangular grids, you cannot
%                    remove points and still use a griddedInterpolant.
%
% generatePlots - boolian. If true, then the decimated data will be
% plotted. Default is false.
%
% Output parameters:
% dcmt - a structure with fields X, Y, and Z. x and y will be the grid
% over which the data is evaluated and z will be the decimated data fit on
% that grid. 
%
% Notes: 
%
% If gridType is 'rectangular', then the Z data will be symmetrically
% padded by on all sides and the X and Y data extended by on all sides.
% This is so that the surface may be used to interpolate the original data
% without producing NaNs on the border.
%
% If the gridType is 'hexagonal', then the output data will be n x 1 column
% arrays where the i'th hexagon will have a center of (X(i),Y(i)) and a
% value of Z(i). Hexagonal grids will not by symmetrically padded like the
% rectangular grids. Hexagonal grids will need to be interpolated with a
% scatteredInterpolant, or you will need to creat a delaunay triangulation
% and then use that to do linear interpolation.
%
% See also CART2HEX HEX2CART


% James Kapaldo
% August 22, 2016
%
% 2016-10-10: Re-named from bg_decimate. Added in density option. Added in
% hexagonal grids. Updated help file.

% Make sure inputs are columns.

if ~iscolumn(x)
    x = x(:);
end
if ~iscolumn(y)
    y = y(:);
end
if ~iscolumn(z)
    z = z(:);
end
z = double(z);

% Get input values
p = inputParser;
p.FunctionName = 'decimateData';

addParameter(p,'binSize',0.25,@(t) numel(t) <= 2 && all(t > 0));
addParameter(p,'reductionMethod','median',@(x) any(strcmpi(x,{'mean','median','mode','percentile','std','mad','density'})) || isa(x,'function_handle'))
addParameter(p,'generatePlots',false,@(t) islogical(t) || (t==0 || t == 1) );
addParameter(p,'percentValue',50,@(t) numel(t) == 1 && t >= 0 && t <= 100 );
addParameter(p,'gridType','rectangular',@(x) any(strcmpi(x,{'rectangular','hexagonal'})));
addParameter(p,'cleanHexagonData',false, @(x) (x==1 || x==0));

parse(p,varargin{:})

voxelSizes = p.Results.binSize;
reductionMethod = lower(p.Results.reductionMethod);
generatePlots = p.Results.generatePlots;
percentValue = p.Results.percentValue;
gridType = lower(p.Results.gridType);
cleanHexagonData = p.Results.cleanHexagonData;

switch gridType
    case 'rectangular'
        if numel(voxelSizes)==1
            voxelSizes(2) = voxelSizes;
        end

        % boundEdges = prctile([x,y],[0,100]);
        boundEdges = [min(x), min(y); max(x), max(y)];
        xEdges = (floor(boundEdges(1,1)/voxelSizes(1))*voxelSizes(1)-voxelSizes(1)) : voxelSizes(1)  : (ceil(boundEdges(2,1)/voxelSizes(1))*voxelSizes(1) + voxelSizes(1)); %this 2 is important, because we need to fit the surface on the bottom so that we have known points all around our actual data, that way when we interpolate the surface we do not get nan's around the edge
        yEdges = (floor(boundEdges(1,2)/voxelSizes(2))*voxelSizes(2)-voxelSizes(2)) : voxelSizes(2)  : (ceil(boundEdges(2,2)/voxelSizes(2))*voxelSizes(2) + voxelSizes(2));

        xCents = xEdges(1:end-1) + diff(xEdges)/2;
        yCents = yEdges(1:end-1) + diff(yEdges)/2;
        [X,Y] = ndgrid(xCents,yCents);

        bin(:,2) = discretize(y,yEdges);
        bin(:,1) = discretize(x,xEdges);
        
    case 'hexagonal'
        if numel(voxelSizes)==1
            voxelSizes(2) = 1;
        end
        
        hx = cart2hex([x,y],voxelSizes(1),voxelSizes(2));

        xEdges = min(hx(:,1))-0.5:max(hx(:,1))+0.5;
        yEdges = min(hx(:,2))-0.5:max(hx(:,2))+0.5;
        xCents = xEdges(1:end-1) + diff(xEdges)/2;
        yCents = yEdges(1:end-1) + diff(yEdges)/2;

        [X,Y] = ndgrid(xCents,yCents);

        bin(:,2) = discretize(hx(:,2),yEdges);
        bin(:,1) = discretize(hx(:,1),xEdges);
end


% n = accumarray(bin(all(bin>0,2),:),1,size(X),[],0);

if ischar(reductionMethod)
    switch reductionMethod
        case 'median'
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@median,nan);
        case 'mode'
            MAD = mad(z,1)/50;
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@(x) mode(round(x/MAD)*MAD),nan);
        case 'mean'
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@mean,nan);
        case 'percentile'
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@(x) prctile(x,percentValue),nan);
        case 'std'
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@std,nan);
        case 'mad'
            Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),@(x) mad(x,1),nan);
        case 'density'
            Z = accumarray(bin(all(bin>0,2),:),1,size(X),[],0);
    end
else
    Z = accumarray(bin(all(bin>0,2),:),z(all(bin>0,2)),size(X),reductionMethod,nan);
end

switch gridType
    case 'rectangular'
        Z(:,1) = Z(:,2);
        Z(1,:) = Z(2,:);

        Z(:,end) = Z(:,end-1);
        Z(end,:) = Z(end-1,:);
    case 'hexagonal'
        
        Z = Z(:);
        sizeX = size(X);
        [X,Y] = hex2cart([X(:),Y(:)],voxelSizes(1),voxelSizes(2));
        
        
        if cleanHexagonData
            if strcmp(reductionMethod,'density')
                toRemove = ~Z;
            else
                n = accumarray(bin(all(bin>0,2),:),1,sizeX,[],0);
                toRemove = ~n(:);
            end
            X(toRemove) = [];
            Y(toRemove) = [];
            Z(toRemove) = [];
        end

end

dcmt.X = X;
dcmt.Y = Y;
dcmt.Z = Z;

if generatePlots
    switch gridType
        case 'rectangular'
            figure
            surface(X,Y,Z)
        case 'hexagonal'
            hexplot([X,Y],voxelSizes(1),voxelSizes(2),'colorData',Z,'sizeData','sameAsColor','colorScale','linear','sizeScale','none','maxHexSize',1);
    end
    
    if ischar(reductionMethod)
        if strcmp(reductionMethod,'percentile')
            title([reductionMethod ' ' sprintf('%0.1f',percentValue) '%'])
        else
            title(reductionMethod)
        end
    else
        title(char(reductionMethod))
    end
    axis tight
    switch gridType
        case 'rectangular'
            colorbar
        case 'hexagonal'
            c = colorbar;
            c.Label.String = ['log( ' char(reductionMethod) ' )'];
            c.FontWeight = 'bold';
    end
    drawnow;
    goDark(gcf)
end

end