function G = interpolateGradient(V,p)
% INTERPOLATEGRADIENT Compute the gradient of V and evaluate it at p.
%
% This function is made to conserve space so that D arrays of size V, where
% D is the dimension, do not need to be saved for the derivative of V along
% each dimension. 
%
% This function is not fast when the number of points, size(p,1), is large.

% James Kapaldo

sz = size(V);
D = length(sz);
N = size(p,1);

G = zeros(N,D);

pl = floor(p);
ph = ceil(p);

pl = min(sz-1,max(1,pl-1));
ph = max(2,min(sz,ph+1));

inds = cell(1,D);
dV = cell(1,D);

for n = 1:N
    
    for i = 1:D
        inds{i} = pl(n,i):ph(n,i);
    end

    Vs = V(inds{:});
    
    [dV{:}] = gradient(Vs);
    dV([1,2]) = dV([2,1]);
    
    for i = 1:D   
        F = griddedInterpolant(inds,dV{i});
        G(n,i) = F(p(n,:));
    end

end



% dV = cell(1,D);
% [dV{:}] = gradient(V);
% dV([1,2]) = dV([2,1]);
% 
% dx = cell(1,D);
% for i = 1:D
%     dx{i} = (1:sz(i));
% end
% 
% for i = D:-1:1
%     F = griddedInterpolant(dx,dV{i});
%     G(:,i) = F(p);
% end


end