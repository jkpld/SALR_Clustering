function [n, cents, sz] = prepare_ND_Example_Data(X, nbins, density_threshold, smooth_data, maintain_aspectRatio)

if nargin < 4
    smooth_data = false;
end
if nargin < 5
    maintain_aspectRatio = false;
end
D = size(X,2); % D dimensions


% Do an initial discretization of the data -------------------------------
boundEdges = [min(X); max(X)];
[n, cents, sz] = discretizeData(X, boundEdges, nbins, 0, 0);


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
[n, cents, sz] = discretizeData(X, boundEdges, nbins, smooth_data, maintain_aspectRatio);

n(n<0) = 0;


end


function [n, cents, sz] = discretizeData(X, boundEdges, nbins, smooth_data,useMean)

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

sz = uint8(cellfun(@length, cents));

if smooth_data 
    n = zeros(sz,'single');
    h = ndGaussianFilter(D,0.5,3);

    bin = bin(all(bin>0,2),:);

    fprintf('Smoothing data ...\n')
    counter = 1;
    for i = -1:1
        for j = -1:1
            fprintf('%d/9 ', counter)
            for k = -1:1
                for u = -1:1
                    for v = -1:1
                        
                        tmp = min(bin + uint8([i,j,k,u,v]),sz);
                        n = n + h(i+2,j+2,k+2,u+2,v+2) * accumarray(tmp,1,sz,[],0);
                        fprintf('.')
                    end
                end
            end
            fprintf('\n')
            counter = counter + 1;
        end
    end
else
    n = accumarray(bin(all(bin>0,2),:),1,sz,[],0);
end

sz = double(sz);

end