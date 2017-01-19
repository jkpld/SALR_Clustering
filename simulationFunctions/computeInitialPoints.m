function [r0,Info] = computeInitialPoints(method,BW,options)
% COMPUTEINITIALPOINTS  Generate a set of initial particle locations for
% modeling.
%
% [r0,Info] = computeInitialPoints(method, BW, options)
%
% Create a set of initial positions r0, using METHOD, from a binary mask BW
% or a set of initial possible positions (see below). OPTIONS should be a
% structure that gives additional information. Each option should be a
% field with the required value. There is one required option for all
% methods, and some methods have additional required options. 
%
% All initial points returned will be located inside of the binary mask.
%
% Input parameters:
% method - The method used to compute the initial points (see below).
% BW - A binary mask of a single object
% options - A structure with additional parameters (see below).
%
% Required option:
% rs - Wigner-Seitz radius giving the effective radius of a particle. This
%      will be used to calculate grid sizes and/or the number of
%      particles.
%
% Available methods:
% 'random' - Select N random locations from the binary mask where N is the
%            area of the mask divided by the effective area of a particle
%            (pi*rs^2).
%
% 'uniform' - Generate a hexagonal grid of points with a lattice constant
%             of 2*rs. Use all lattice sites that are inside the binary
%             mask as the set of initial positions.
%
% 'uniformRandom' - Overlay a hexagonal grid over the binary mask with a
%                   lattice constant of 2*rs. In each hexagon, choose one
%                   random point from the mask.
%
% 'r0set_random' - From the set of possible initial positions, choose N
%                  random points, where N is the area of the mask divided
%                  by the effective particle area (pi*rs^2).
%                       
% 'r0set_uniformRandom' - Generate a hexagonal grid with lattice constant
%                         2*rs. In each hexagon, choose one random position
%                         from the set of possible initial positions
%
% Additional requierd options for each method:
% There are no additional required options for 'random' or 'uniform'. For
% both 'r0set_random' and 'r0set_uniformRandom' there is one
% additional required options
%
% r0set - An array of possible initial positions.
%
% Output parameters:
% r0 - Nx2 array with the initial N particle locations.
% Info - A structure array giving information that may be off help when
%        debuging. Fields:
%           N - the number of particles returned (number of rows in r0).
%
%           givenMethodFailed - boolian. If true, then method requested
%           failed to produce a single particle position. (See notes below
%           for more information about this.)
%
%           r0setHexData - A structure with fields X, Y, Z with the
%           data returned from using decimateData on the curvature centers
%           returned by getCurvatureAndShapeMarkers. This field is only
%           returned if the method is 'r0set_uniformRandom'.
%
% Notes:
%
% * This function is run for a single object, not an image with many
% objects.
%
% * If the methods 'uniform', 'r0set_random', or
% 'r0set_uniformRandom' fail to produce a single particle location, then
% the 'random' method will be used.
%
% * If the binary mask area is smaller than the effective particle area,
% then one point will still be output.
%
% See also DECIMATEDATA COMPUTEBOUNDARYINFORMATION

% JamesKapaldo
% 2016-10-10


% There must be an option giving the Wigner-Seitz radius of each particle
if ~isfield(options,'rs')
    error('getInitialPoints:missingInput','You must have a field named "rs" in OPTIONS that gives the Wigner-Seitz radius of each particle (if the particle where a sphere/circle this would be its radius).')
end

% Get rs.
rs = options.rs;

% Compute effective area of particle.
A = pi*rs^2;

% Get size of image.
imSize = size(BW);

% Initialize givenMethodFailed if Info is requested.
if nargout > 1
    Info.givenMethodFailed = false;
end

% Get initial information for the methods and test to make sure a valid
% method was given
if ~any(strcmp(method,{'random','uniform','uniformRandom','r0set_random','r0set_uniformRandom'}))
    error('getInitialPoints:unknownMethod','Unknown method. Allowed methods are ''random'', ''uniform'', ''uniformRandom'', ''r0set_random'', ''r0set_uniformRandom''.')
end

if any(strcmp(method,{'random','r0set_random'}))
    % The number of particles N is given by the area of the mask
    % divided by the area of a particle

    N = max(round(sum(BW(:))/A),1);
end

if any(strcmp(method,{'uniform','r0set_uniformRandom','uniformRandom'}))
    % Hexagonal grid lattice constant
    a = 2*rs;
end

if any(strcmp(method,{'r0set_random','r0set_uniformRandom'}))
    % There must be an option giving the set of initial possible positions
    if ~isfield(options,'r0set')
        error('getInitialPoints:missingInput','You must have a field named "r0set" in OPTIONS that curvature centers.')
    end

    markers = options.r0set;
    markers(markers(:,1) < 1 | markers(:,1) > imSize(1) | markers(:,2) < 1 | markers(:,2) > imSize(2),:) = [];
    % Remove markers outside the mask.
    markers_lin = markers(:,1) + (markers(:,2)-1)*imSize(1);
    markers(~BW(markers_lin),:) = [];
    
end


% Create the initial points
switch method
    case 'random'
        % Select N random points from the BW mask.
        
        % Use the x-y locations of the eroded mask as the start points. The
        % mask is first eroded because we do not want any start point close
        % to the boundary (where the particle would start with high
        % potential energy).
        [i,j] = find(BW);
        validInds = [i,j];
        
        if size(validInds,1) < N
            warning('getInitialPoints:N_to_large','The number of possible starting locations is smaller than the requeted number of points. Using all possible starting locations.')
            N = size(validInds,1);
        end
        
        % Choose a random N start points.
        r0 = validInds(randperm(size(validInds,1),N),:);
        
    case 'uniform'
        % Create a hexagonal grid of points with lattice constant of twice
        % the Wigner-Seitz radius and take all points inside of the
        % (eroded) mask.
        
        % Create hexagonal grid.
        x = a/2 : a : imSize(2)+a;
        y = a/2 : sqrt(3)*a/2 : imSize(1)+a;
        [X,Y] = ndgrid(x,y);

        X(:,2:2:end) = X(:,2:2:end) + a/2;
        
        % Create start points.
        r0 = round([Y(:),X(:)]);
        r0( r0(:,1) > imSize(1) | r0(:,2) > imSize(2) ,:) = [];
        
        % Get the indices of the starting points that are inside of the
        % mask.
        good = BW( r0(:,1) + (r0(:,2)-1)*imSize(1) );
        
        % Remove all starting points not inside the mask.
        r0 = r0(good,:);
        
        % Ensure that if there were not centers inside the mask we still
        % get some results
        if isempty(r0)
            r0 = computeInitialPoints('random',BW,B,options);
            if nargout > 1
                Info.givenMethodFailed = true;
            end
        end
    case 'uniformRandom'
        % Overlay a hexagonal grid with lattice constant 2*rs and then
        % select one pixel in the mask from each hexagon.
        
        % Put this in a try-catch structure to handle the case when the
        % region is very small in one direction. This case would give an
        % error when we are creating a hexagonal grid with having only one
        % bin edges (when two bin edges are required)
        try
            [i,j] = find(BW);
            markers = [i,j];
            mu = mean(markers,1) - rs/2;
            markers = markers-mu;

            dcmt = decimateData(markers(:,1),markers(:,2),1:size(markers,1),...
                'reductionMethod',@(x) x(randi(numel(x),1)),...
                'gridType','hexagonal',...
                'binSize',a,'cleanHexagonData',1);

            r0_inds = dcmt.Z(~isnan(dcmt.Z));
            r0 = markers(r0_inds,:) + mu;

            % Ensure that if there were not centers inside the mask we still
            % get some results
            if isempty(r0)
                r0 = computeInitialPoints('random',BW,B,options);
                if nargout > 1
                    Info.givenMethodFailed = true;
                end
            end

            if nargout > 1
                dcmt.X = dcmt.X + mu(1);
                dcmt.Y = dcmt.Y + mu(2);
                Info.r0setHexData = dcmt;
            end
        catch ME
            if strcmp(ME.identifier,'MATLAB:discretize:EmptyOrScalarEdges')
                r0 = computeInitialPoints('random',BW,B,options);
                if nargout > 1
                    Info.givenMethodFailed = true;
                end
            else
                rethrow(ME)
            end
        end
    case 'r0set_random'
        % From the curvature, get the radius and the center of the circle.
        % Use a random N of the centers as the initial points.
                
        markers = unique(markers,'rows');
        
        if size(markers,1) < N
            warning('getInitialPoints:N_to_large','The number of possible starting locations is smaller than the requeted number of points. Using all possible starting locations.')
            N = size(markers,1);
        end
        
        r0 = markers(randperm(size(markers,1),N),:);
        
        % Ensure that if there were not centers inside the mask we still
        % get some results
        if isempty(r0)
            r0 = computeInitialPoints('random',BW,B,options);
            if nargout > 1
                Info.givenMethodFailed = true;
            end
        end
        
    case 'r0set_uniformRandom'
        % From the curvature, get the radius and the center of the circle.
        % Overlay a hexagonal grid with lattice constant 2*rs and then
        % select one curvature center from each hexagon.
        
        % Put this in a try-catch structure to handle the case when the
        % region is very small in one direction. This case would give an
        % error when we are creating a hexagonal grid with having only one
        % bin edges (when two bin edges are required)
        try
            [i,j] = find(BW);
            mu = [mean(i),mean(j)] - rs/2;
            markers = markers-mu;

            dcmt = decimateData(markers(:,1),markers(:,2),1:size(markers,1),...
                'reductionMethod',@(x) x(randi(numel(x),1)),...
                'gridType','hexagonal',...
                'binSize',a,'cleanHexagonData',1);

            r0_inds = dcmt.Z(~isnan(dcmt.Z));
            r0 = markers(r0_inds,:) + mu;

            % Ensure that if there were not centers inside the mask we still
            % get some results
            if isempty(r0)
                r0 = computeInitialPoints('random',BW,B,options);
                if nargout > 1
                    Info.givenMethodFailed = true;
                end
            end

            if nargout > 1
                dcmt.X = dcmt.X + mu(1);
                dcmt.Y = dcmt.Y + mu(2);
                Info.r0setHexData = dcmt;
            end
        catch ME
            if strcmp(ME.identifier,'MATLAB:discretize:EmptyOrScalarEdges')
                r0 = computeInitialPoints('random',BW,B,options);
                if nargout > 1
                    Info.givenMethodFailed = true;
                end
            else
                rethrow(ME)
            end
        end
    otherwise
        error('getInitialPoints:unknownMethod','Unknown method. Allowed methods are ''random'', ''uniform'', ''r0set_random'', ''r0set_uniformRandom''. (In switch case.)')
end % switch

if nargout > 1
    Info.N = size(r0,1);
end

end
% computeInitialPoints : changeLog
% 2016-11-2 : added in uniformRandom method
% 2017-01-20 : renamed options and removed the eroding and potential limit
% of the mask