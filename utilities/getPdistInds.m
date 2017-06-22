function inds = getPdistInds(N)
% GETPDISTINDS  Get the indices of the pairs formed by pdist.

% James Kapaldo
% 2016-09-27

if N < 350
    [X,Y] = meshgrid(1:N,1:N);
    mask = tril(true(N),-1);
    
    inds = [X(mask),Y(mask)];
else
    inds = zeros(N*(N-1)/2,2);
    
    tmpI = (1:N)';
    tmp1 = ones(N,1);
    counter = 1;
    for i = 1:N-1
        
        inds(counter:counter+N-i-1,:) = [i*tmp1(i+1:N), tmpI(i+1:N)];
        counter = counter + N-i;
        
    end
end

end
    
