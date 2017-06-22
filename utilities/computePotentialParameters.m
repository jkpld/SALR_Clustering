function computePotentialParameters(depths,minimum_Locations,extents,file)
% COMPUTEPOTENTIALPARAMETERS  Solve for the potential parameters that give
% all combinations of the DEPTHS, MINIMUM_LOCATIONS, and EXTENTS given.
% Save the paramters in FILE

% James Kapaldo
% 2016-10-11
% 2016-10-28 : added in the input parameters

% depths = -0.1:-0.1:-2;
% centers = 1.5:0.5:4;
% extents = 8:1:20;

centers = minimum_Locations;

potentialParameters = zeros(numel(depths),numel(centers),numel(extents),3);

options = optimset('display','none','TolX',1e-7,'TolFun',1e-7,'MaxFunEvals',500);

for i1 = 1:numel(depths)
    for i2 = 1:numel(centers)
        for i3 = 1:numel(extents)
            r = centers(i2)/2:0.01:(extents(i3) + centers(i2));
            f = @(x) potentialParametersFun(x,[depths(i1),centers(i2),extents(i3)],r);
            potentialParameters(i1,i2,i3,:) = fminsearch(f,[-depths(i1),centers(i2),(extents(i3)-centers(i2))/3],options); % 
        end
    end
    fprintf('%d/%d...\n',i1,numel(depths))
end

PotentialParameters.depth = depths;
PotentialParameters.center = centers;
PotentialParameters.extent = extents;
PotentialParameters.parameters = potentialParameters; %#ok<STRNU>

% potParams = PotentialParameters.parameters;

save(file,'-struct','PotentialParameters')

end

function f = potentialParametersFun(x,y,r)

Vint = 1./(r+0.2) - x(1)*exp(-(r-x(2)).^2/(2*x(3)^2));

[minV,minVind] = min(Vint);

[~,extentInd] = max(Vint(minVind:end));

f = sum((([minV,r(minVind),r(extentInd+minVind-1)] - y)./y).^2,2);

end