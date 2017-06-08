function dV = confiningForce(V, problem_scales, options)
% CONINFINGFORCE : Return function handle that takes in particle positions
% and returns the force from the confining potential.
%
% dV = confiningForce(V, problem_scales, options)
%
% Input parameters,
% V : confining potential
% problem_scales : structure with the two fields,
%   grid_spacing : 1xD array giving the spacing of the confining
%       potential's grid. D is the dimension of the data.
%   grid_to_solver : 1xD array giving the scaling factors to convert from
%       the grid to the solver space.
% options : instance of class seedPointOptions
%
% Output parameters,
% dV : function handle that takes in an array of positions (NxD, where D is
%   the dimension) and outputs the gradients at those positions.
%   Example. p = dV(r), p will have the same size as r.
%
%   Note, the positions passed to dV() should be in the solver space. This
%   function will automatically convert the positions back into grid space
%   to evaluate the gradients.

% James Kapaldo

grid_spacing = problem_scales.grid_spacing;
grid_to_solver = problem_scales.grid_to_solver;

MAX_MEMORY = options.Maximum_Memory;
MAX_FORCE = options.Max_Potential_Force;

sz = size(V);
D = length(sz);
solver_to_grid = 1./grid_to_solver;

if sizeof(class(V))*D*numel(V)/1e9 > MAX_MEMORY
    % If the amount of space taken up by all of the gradients is more than
    % MAXIMUM_MEMRORY GB, the change how the derivatives are computed to
    % save space.
    
    if ~isnan(MAX_FORCE)
        str = ['The memory requirement is above the Maximum_Memory ',...
            'threshold and the Max_Potential_Force has been set. The ',...
            'current implimentation of the memory efficent force ',...
            'calculation does not allow of theforce to be scaled. ',...
            'You will need to impliment this yourself.'];
        
        warning('confiningForce:maxForceSet',str)
    end
    
    dV = @(r) interpolateGradient(V, grid_spacing, solver_to_grid, r);
else
    % Compute gradient of confining potential
    dVi = cell(1,D);
    dtscl = num2cell(grid_spacing([2,1,3:end]));
    [dVi{:}] = gradient(V, dtscl{:});
    dVi([1,2]) = dVi([2,1]);
    % 1 and 2 are flipped above because gradient() always outputs the
    % derivative along the second dimension first.
    
    % Normalize the gradient magnitude if required
    if ~isnan(MAX_FORCE)
        M = zeros(size(V),'like',V);
        
        % Compute gradient magnitude
        for i = 1:D
            M = M + dVi{i}.^2;
        end
        
        p99 = sqrt(prctile(M(V<1),99)); % 99th percentile of gradient magnitude inside the object.
        
        % Scale the gradients
        for i = 1:D
            dVi{i} = MAX_FORCE * dVi{i} / p99;
        end
    end
    
    if D == 2
        for i = D:-1:1
            dVi{i} = @(r) interp2mex(dVi{i}, r(:,2), r(:,1));
            % This mex interpolation function assumes the grid spacing is
            % 1, which is what our spacing is in grid space, and it uses
            % nearest neighbor extrapolation.
        end
        % Note, if the above mex function is not installed, then just use
        % the code below in the "else" case. It will work fine, but will be
        % a bit slower.
    else
        % Create grid vectors
        for i = D:-1:1
            dx{i} = (1:sz(i));
        end
        % Create griddedInterpolants
        for i = D:-1:1
            dVi{i} = griddedInterpolant(dx,dVi{i});
        end
    end
    
    dV = @(r) formGradient(dVi,r.*solver_to_grid);
end

end

function G = formGradient(dVi, r)
D = length(dVi);
G = zeros(size(r,1),D);
for i = 1:D
    G(:,i) = dVi{i}(r);
end
end

function G = interpolateGradient(V, grid_spacing, solver_to_grid, p)
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

p = p .* solver_to_grid;

pl = floor(p);
ph = ceil(p);

pl = min(sz-1,max(1,pl-1));
ph = max(2,min(sz,ph+1));

inds = cell(1,D);
dV = cell(1,D);

dtscl = num2cell(grid_spacing([2,1,3:end]));

for n = 1:N
    
    for i = 1:D
        inds{i} =  pl(n,i):ph(n,i);
    end
    
    Vs = V(inds{:});
    
    [dV{:}] = gradient(Vs, dtscl);
    dV([1,2]) = dV([2,1]);
    
    for i = 1:D
        F = griddedInterpolant(inds,dV{i});
        G(n,i) = F(p(n,:));
    end
    
end

end