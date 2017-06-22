function [B,n,curvature,curvatureCenters] = computeBoundaryInformation(BW,options)
% COMPUTEBOUNDARYINFORMATION  Compute the boundary contours, inward
% pointing normal vectors, curvature, and curvature centers.
%
% [B,n,kappa,curvatureCenters] = computeBoundaryInformation(BW,k)
%
% Take in a binary image BW and compute the boundaries with holes using
% bwboundaries(). Compute the curvature of each boundary contour, compute
% the "curvature centers" (see output parameters below for description),
% and compute the inward pointing boundary normals. Finally, combine outer
% contours and holes into a single array using a row of NaN delimiters to
% seperate individual contours. Any hole with an circumference less than
% the circumference of a circle with area MINHOLEAREA will be removed.
%
% Input parameters:
% BW - binary mask from which the boundaries will be computed
% options - must have the fields
%           Curvature_Smoothing_Size - The size of the smoothing filter
%           used to compute the curvatures. Use integer values. Perhaps try
%           2 too see how it works first.
%
%           Max_Radius - This value will set a threshold for what is
%           considered positive curvature (regions of positive curvature
%           are concave regions of the boundary). The threshold will be set
%           as 1/(2*maximumRadius).
%
%           Use_Parallel - logical. If true, then a parallel loop will be
%           used to calculate the curvature information.
%
% Output parameters:
% B - A cell array where each element (Mx2 array) has a parent (outer) 
%     boundary followed bay any holes with a row of NaN delimiters.
%       Example. If the i'th parent boundary has N holes, then B{i} will
%       have the form:
%             B{i} = [ parentContour; 
%                      NaN, NaN; 
%                      childContour1; 
%                      NaN, NaN;
%                      childContour2;
%                      NaN, NaN;
%                      ...
%                      NaN, NaN;
%                      childContourN ];
%
% n - A cell array where each element (Mx2 array) has the inward pointing
%     boundary normals for the contours. n has the same format as B.
%
% curvature - A cell array where each element (Mx1 array) has the 
%             curvatures for the contours. kappa has the same format as B.
%
% curvatureCenters - A cell array where each element contains the centers
%                    of the circles computed from the curvatures and the
%                    inward pointing boundary normals.
%
% Notes:
% Parent contours will be oriented clockwise, children contours (holes)
% will be oriented counter-clockwise.
%
% See also GETCURVATUREANDSHAPEMARKERS

% James Kapaldo
% 2016-10-10
% 2016-10-28 : added in parallel version. It saves a bit of time, on one
%              image it goes from ~1.9 sec to ~1.5 sec. bwboundaries takes
%              up almost all of the time in this function (1.3 sec for the
%              image used.) Changed input to options structure (class
%              declumpOptions)


% USE_PARALLEL = options.Use_Parallel;
KAPPA_SMOOTHING_SIGMA = options.Curvature_Smoothing_Size;
MAX_RADIUS = options.Max_Radius;

% if USE_PARALLEL
%     [B,n,curvature,curvatureCenters] = computeBoundaryInformation__parallel(BW,KAPPA_SMOOTHING_SIGMA,MAX_RADIUS);
%     return;
% end

useConvexHull = true;

% Get the boundaries and parant-child matrix
[tmpB,~,numObjs,bndryTplgy] = bwboundaries(BW,8);

% Get the image size
imSize = size(BW);

% Initialize cell arrays
B = cell(1,numObjs);
n = cell(1,numObjs);
curvatureCenters = cell(1,numObjs);
curvature = cell(1,numObjs);

% Iterate over each object finding the curatures, normals, and curvature
% centers of the object and its children (if any).
for i = 1:numObjs
    
    % Get parent contour and oriented it clockwise
    Prnt = tmpB{i};
    if ~ispolycw(tmpB{i}(:,1),tmpB{i}(:,2))
        Prnt = flip(Prnt,1);
    end
    
    % Get curvature, tangents, and markers
    [kappaP,~,~,nP,markersP] = getCurvatureAndShapeMarkers(Prnt,imSize,KAPPA_SMOOTHING_SIGMA,MAX_RADIUS,useConvexHull);
    Prnt(end,:) = [];%nan;
    
    % Find any children
    children = bndryTplgy(:,i);
    if any(children)
        
        % If there are children, then iterate over each one computing the
        % same things.
        children = find(children);
        numChildren = numel(children);
        
        child = cell(1,numChildren);
        kappa_child = cell(1,numChildren);
        n_child = cell(1,numChildren);
        markers_child = cell(1,numChildren);
        
        for j = 1:numChildren

            child{j} = tmpB{children(j)};
            if ispolycw(child{j}(:,1),child{j}(:,2))
                child{j} = flip(child{j});
            end

            [kappa_child{j},~,~,n_child{j},markers_child{j}] = getCurvatureAndShapeMarkers(child{j},imSize,KAPPA_SMOOTHING_SIGMA,MAX_RADIUS,useConvexHull);
            child{j}(end,:) = nan;

        end
        
        % Put the results of the children together with the parent
        B{i} = [Prnt; nan, nan; cat(1,child{:})];
        curvature{i} = [kappaP; cat(1,kappa_child{:})];
        n{i} = [nP; cat(1,n_child{:})];
        curvatureCenters{i} = [markersP; cat(1,markers_child{:})];
        
    else
        B{i} = [Prnt; nan, nan];
        curvature{i} = kappaP;
        n{i} = nP;
        curvatureCenters{i} = markersP;
    end
    
    % Remove the row of NaN's at the end.
    B{i}(end,:) = [];
    curvature{i}(end,:) = [];
    n{i}(end,:) = [];
    
    % Note getCurvatureAndShapeMarkers returns the boundary tangents, so we
    % need to convert them to the inward pointing boundary normals.
    n{i} = -[-n{i}(:,2), n{i}(:,1)]; % Mx2
    n{i} = n{i} ./ sqrt(sum(n{i}.^2,2)); % Mx2
    
end % for

end


function [B,n,curvature,curvatureCenters] = computeBoundaryInformation__parallel(BW,k,maxR)

KAPPA_SMOOTHING_SIGMA = k;
MAX_RADIUS = maxR;
useConvexHull = true;

% Get the boundaries and parant-child matrix
[B_full,~,numObjs,bndryTplgy] = bwboundaries(BW,8);
N = numel(B_full);

% Get the image size
imSize = size(BW);

% Oriente the boundaries, get the boundary curvatures, and get the centers
% of curvature
curvatureCenters_full = cell(1,N);
curvature_full = cell(1,N);
normals_full = cell(1,N);
% B_full = B_full;

parfor contour = 1:N
    if contour <= numObjs
        if ~ispolycw(B_full{contour}(:,1),B_full{contour}(:,2))
            B_full{contour} = flip(B_full{contour},1);
        end
    else
        if ispolycw(B_full{contour}(:,1),B_full{contour}(:,2))
            B_full{contour} = flip(B_full{contour},1);
        end
    end
    % Get curvature, tangents, and markers
    [curvature_full{contour},~,~,normals_full{contour},curvatureCenters_full{contour}] = getCurvatureAndShapeMarkers(B_full{contour},imSize,KAPPA_SMOOTHING_SIGMA,MAX_RADIUS,useConvexHull);

    B_full{contour}(end,:) = [nan,nan]; 
    
    % Note getCurvatureAndShapeMarkers returns the boundary tangents, so we
    % need to convert them to the inward pointing boundary normals.
    normals_full{contour} = -[-normals_full{contour}(:,2), normals_full{contour}(:,1)]; % Mx2
    normals_full{contour} = normals_full{contour} ./ sqrt(sum(normals_full{contour}.^2,2)); % Mx2
end

% Initialize cell arrays
B = cell(1,numObjs);
n = cell(1,numObjs);
curvatureCenters = cell(1,numObjs);
curvature = cell(1,numObjs);

% Iterate over each object finding the curatures, normals, and curvature
% centers of the object and its children (if any).
for i = 1:numObjs
    % Find any children
    children = bndryTplgy(:,i);
        
    % If there are children, then iterate over each one computing the
    % same things.
    B{i} = [B_full{i}; cat(1,B_full{children})];
    curvature{i} = [curvature_full{i}; cat(1,curvature_full{children})];
    n{i} = [normals_full{i}; cat(1,normals_full{children})];
    curvatureCenters{i} = [curvatureCenters_full{i}; cat(1,curvatureCenters_full{children})];
    
    % Remove the row of NaN's at the end.
    B{i}(end,:) = [];
    curvature{i}(end,:) = [];
    n{i}(end,:) = [];    
end % for


end