function inds = getPdistInds(N)
% GETPDISTINDS Get the indices of the pairs formed by pdist().
%
% inds = getPdistInds(N)
%
% Input parameters:
% N : The number of objects.
%
% Output parameters:
% inds : N*(N-1)/2 x 2 array where the ith row gives the indices of the
%   objects for the ith distance returned by pdist().

% Ex.
% x = rand(4,2);
% d = pdist(x);
% inds = getPdistInds(4);
% fprintf('The distance between point %d and %d is %f\n',inds(3,1),inds(3,2),d(3))

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
    


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
