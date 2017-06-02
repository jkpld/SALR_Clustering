function [kappa,CoC_linIdx,negKappaInds,d_bndry,CoC,kappa_s] = getCurvatureAndCenters(BWperim,imSize,kappaSmoothingSigma,maxRadius,useConvexHull)
% GETCURVATUREANDCENTERS  Compute boundary curvature and centers.
%
% [kappa,CoC_linIdx,negKappaInds,d_bndry,CoC,kappa_s] = getCurvatureAndCenters(BWperim,imSize,kappaSmoothingSigma,maxRadius,useConvexHull)

% James Kapaldo


% kappaSmoothingSigma = 4;
negativeMADKappaMultiplier = 2;

% kappaSmoothingSigma is the std of the gaussian smoothing filter applied
% to the boundary before the derivatives are calculted. The derivatives
% will be calculated using the derivative of the gaussian with this same
% std.

% negativeMADKappaMultiplier sets the threshold for determining if a point
% should be considered as a curvature maximum. The threshold will be
% calculated using this through:
%  threshold = negativeMADKappaMultiplier * mad(kappa,1) + median(kappa)

if nargin < 5
    useConvexHull = true;
else
    useConvexHull = logical(useConvexHull);
end


minR = 1;
minK = 1/minR;

POSITIVE_CURVATURE_THRESHOLD = 1/(2*maxRadius);

bndry = BWperim;
cntr = mean(bndry,1);

bndry = bsxfun(@minus,bndry,cntr);

% Create filters for getting derivatives
filtSize = round(7*kappaSmoothingSigma);
x = -floor(filtSize/2):ceil(filtSize/2);
G = exp(-(x).^2/(2*kappaSmoothingSigma^2))/(kappaSmoothingSigma*sqrt(2*pi));
dG = -(x) .* G / kappaSmoothingSigma^2;
ddG = - G / kappaSmoothingSigma^2 + (x).^2 .* G / kappaSmoothingSigma^4;

G = G(:);
dG = dG(:);
ddG = ddG(:);


% Get boundary curvature
bndrys = imfilter(bndry,G,'circular','conv','same');
d_bndry = imfilter(bndrys,dG,'circular','conv','same');
dd_bndry = imfilter(bndrys,ddG,'circular','conv','same');

bndry(end,:) = [];
d_bndry(end,:) = [];
dd_bndry(end,:) = [];

kappa = (d_bndry(:,1).*dd_bndry(:,2) - d_bndry(:,2).*dd_bndry(:,1))./(sum(d_bndry.^2,2).^(3/2));

kappa(kappa < -minK) = -minK; % Changed from "kappa(kappa < -minK) = minK;" on 2016-10-01
kappa(isnan(kappa)) = 0;

% Get boundary normals
n = [-d_bndry(:,2), d_bndry(:,1)];
n = bsxfun(@rdivide,n,sqrt(sum(n.^2,2)));

bndry = round(bsxfun(@plus,bndry,cntr));

% Smooth the negative curvatures out using a convex hull. -- this will make
% the curvature along the long side of an ellipse the same as the curvature
% on the short side, and it will allow of the the shape markers to better
% separate two touching ellipses.

CoC = [];
shiftInd = find(kappa>POSITIVE_CURVATURE_THRESHOLD,1);
kappa_s = kappa;

if ~isempty(shiftInd)

    if useConvexHull
        kappa = circshift(kappa,-shiftInd+1);
        kappa_s = circshift(kappa_s,-shiftInd+1);
        nInd = kappa < POSITIVE_CURVATURE_THRESHOLD;
        % Remove single positive (negative) points surrounded by negative
        % (positive) points : remove isolated positive or negative points
        nInd = sum([nInd,circshift(nInd,-1,1),circshift(nInd,1,1)],2)>1;

        bnds = [0;find(abs(diff(nInd)));numel(kappa)];
        ind = (1:numel(kappa))';

        for i = 2:2:numel(bnds)-1
            % only replace the negative curvatures with the smoothed convex
            % hull. can also smooth the positive curvatures by replacing with
            % 1:numel(bnds)-1, but then finding the maximum curvatures for
            % delunay triangulation does not work as well.
            tmpind = bnds(i)+1:bnds(i+1);
            if numel(tmpind) < 2
                kappa_s(tmpind) = kappa(tmpind);
            else
                if mod(i,2)==1

                    tmpx = [ind(bnds(i)+1); ind(tmpind); ind(bnds(i+1))];
                    tmpy = [-2;kappa(tmpind);-2];
                    K = convhull(tmpx,tmpy);
                    K(tmpy(K)==-2) = [];

                else
                    tmpx = [ind(bnds(i)+1); ind(tmpind); ind(bnds(i+1))];
                    tmpy = [2;kappa(tmpind);2];
                    K = convhull(tmpx,tmpy);
                    K(tmpy(K)==2) = [];

                end

                kappa_s(bnds(i)+1:bnds(i+1)) = nakeinterp1(tmpx(K),tmpy(K),ind(tmpind));

            end
        end

        kappa_s = circshift(kappa_s,shiftInd-1);
        kappa = circshift(kappa,shiftInd-1);

    end

    R_s = 1./kappa_s;
    CoC = bndry + bsxfun(@times, n,R_s);

end

% Get the indices of the positive curvature peaks for use with delunay
% triagulation.
% negKappaInds = find((kappa_s == imdilate(kappa_s,ones(5,1))) & kappa_s > max((negativeMADKappaMultiplier*mad(kappa_s,1)+median(kappa_s)),0));
negKappaInds = [];

if isempty(CoC)
    CoC = bsxfun(@plus, cntr, [-1 0; 1 0; 0 -1; 0 1; 0 0]);
end

CoC = round(CoC);

% convert the centers to linear indices of the full image.
CoC_linIdx = CoC(:,1) + (CoC(:,2)-1)*imSize(1);

% set the output kappa, complete the contour and add a row of nan.

kappa = [kappa; nan];
kappa_s = [kappa_s; nan];
d_bndry = [d_bndry; nan, nan];



end
