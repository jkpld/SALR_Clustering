function [n, cents, sz] = binData(X, nbins, density_threshold, smooth_data, maintain_aspectRatio)

if nargin < 4
    smooth_data = false;
end
if nargin < 5
    maintain_aspectRatio = false;
end

D = size(X,2); % D dimensions

% Do an initial discretization of the data -------------------------------
boundEdges = [min(X); max(X)];
[n, cents, sz] = discretizeData(X, boundEdges, nbins, 0);

% Find range where data is above density threshold ----------------------
idx = find(n >= density_threshold);

locs = cell(1,D);
[locs{:}] = ind2sub(sz,idx);
idx = cat(2,locs{:});

minIdx = max(1,min(idx)-1);
maxIdx = min(sz,max(idx)+1);

for i = 1:D
    boundEdges(:,i) = [cents{i}(minIdx(i)), cents{i}(maxIdx(i))];
end

% Re-discretize the data with the same number of bins across the region 
% above the density threshold. ------------------------------------------
[n, cents, sz] = discretizeData(X, boundEdges, nbins, maintain_aspectRatio);

if smooth_data
    n = smooth_nd(n, 0.5, 3, 1);
end

n(n<0) = 0;

end

function [n, cents, sz] = discretizeData(X, boundEdges, nbins, useMean)

D = size(X,2);
vS = diff(boundEdges,1)/nbins;
if useMean % Maintain data aspect ratio
    vS = mean(vS)*ones(1,D);
end

cents = cell(1,D);
bin = zeros(size(X,1),D,'uint8');
for i = D:-1:1
    edges = (floor(boundEdges(1,i)/vS(i))*vS(i)-vS(i)) : vS(i)  : (ceil(boundEdges(2,i)/vS(i))*vS(i) + vS(i)); 
    bin(:,i) = discretize(X(:,i),edges);    
    cents{i} = edges(1:end-1) + vS(i)/2;
end

sz = cellfun(@length, cents);
n = accumarray(bin(all(bin>0,2),:),1,sz,[],0);
end

function n = smooth_nd(dat, sigma, hsize, verbose)

    if nargin < 4
        verbose = false;
    end

    % Get the indices of the non-zero data
    sz = int16(size(dat));
    D = length(sz);
    
    dat_lin_idx = find(dat);
    [dat_idx{1:D}] = ind2sub(sz,dat_lin_idx);
    dat_idx = int16(cat(2,dat_idx{:}));   
    
    % Create guassian filter and index arrays
    h = ndGaussianFilter(D,sigma,hsize);
    r = (hsize-1)/2; % filter radius

    % Get the index arrays for each element in h
    idx = int16(perm_with_rep(-r:r,D));
    for ii = D:-1:1
        lin_idx{ii} = idx(:,ii) + r + 1;
    end
    lin_idx = sub2ind(size(h),lin_idx{:});

    
    % Initialize the output
    n = zeros(sz,'single');
    
    if verbose
        fprintf('\nSmoothing data...   0%%\n')
    end
    
    % Iterate over each element in the smoothing filter
    for ii = 1:hsize^D
        % Offset the data indices by filter index (filter indices are
        % centered on zero)
        tmp = min(max(dat_idx + idx(ii,:),1),sz);        
        
        % Compute the new linear indices
        tmp_lin = sub2ind_fast(sz, tmp);

        % Add the data, scaled by the filter, to the offset indices of the
        % output
        n(tmp_lin) = n(tmp_lin) + h(lin_idx(ii)) * dat(dat_lin_idx);

        if verbose
            fprintf('\b\b\b\b\b%3.0f%%\n',100*ii/hsize^D)
        end
    end
    
    if verbose
        fprintf('\b ...Finished!\n')
    end

end

function out = perm_with_rep(V,N)

% https://www.mathworks.com/matlabcentral/fileexchange/7147-permn-v--n--k-

nV = numel(V);
[Y{N:-1:1}] = ndgrid(1:nV);
I = reshape(cat(N+1,Y{:}),[],N);
out = V(I);

end

function lin_inds = sub2ind_fast(sz, inds)

sz = double(sz);
sz = cumprod(sz);

lin_inds = double(inds(:,1));
for i = 2:length(sz)
    lin_inds = lin_inds + (double(inds(:,i))-1)*sz(i-1);
end

% The above is much faster than the code below.
% sz = circshift(sz,1);
% sz(1) = 1;
% lin_inds = (double(inds) + [0,-1,-1,-1,-1])*sz';

end