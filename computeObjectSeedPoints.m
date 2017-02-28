function [seedPoints, Info] = computeObjectSeedPoints(BW, M, r0set, useCentroid, options, objNumber, errorCount)

% Info.Message : 0, everything is fine
%                1, object is convex or too small (the center of the object will be the seed point)
%                2, less than 2 initial particles (the center of the object will be the seed point)
%                3, there was twice an error, object will be skipped

% See also MODELPARTICLEDYNAMICS EXTRACTCLUSTERCENTERS COMPUTEINITIALPOINTS
% ADD_EXTERIOR_CONFINING_POTENTIAL PROCESSOBJECTS

% M : base potential modifier - should be 1 outside of nuclei.


Info = [];

% Get the error count
if nargin == 6
    errorCount = 0;
end

DEBUG = options.Debug;

if useCentroid
    [seedPoints, Info] = computeCentroid(BW,DEBUG,1);%'convex_or_tooSmall');
    return
end % if

try

    PAD_SIZE = options.Potential_Padding_Size;
    R0_MAXV = options.Maximum_Initial_Potential;
    POINT_SELECTION_METHOD = options.Point_Selection_Method;

    % Create full confining potential ---------------------------------------
    [V, scaleFactor, overlapFactor] = create_scaleInvar_confining_potential(BW, options);
    options.Scale_Factor = scaleFactor;

    % At this point you could use the overlapFactor, the 3rd output of
    % create_scaleInvar_confining_potential, to modify some parameters,
    % for example, the Potential_Extent.

    if ~isempty(M)
        V = V .* M; % Add in any modifier
    end

    % Get initial particle locations ----------------------------------------
    initPointOptions.rs = options.Wigner_Seitz_Radius;
    initPointOptions.r0set = r0set;

    % Don't let the particles start with too high a potential energy.
    allowed_r0_mask = ( V(PAD_SIZE+1:end-PAD_SIZE, PAD_SIZE+1:end-PAD_SIZE) < R0_MAXV ); % & imerode(BW,strel('disk',3));

    if ~any(allowed_r0_mask(:))
        if ~any(BW)
            error('computeObjectSeedPoints:zeroMask','The object mask is all 0''s! I don''t know how this could happen.')
        else
            warning('computeObjectSeedPoints:noValidPositions','There are no valid positions to put a particle. Try increasing the Maximum_Initial_Potential.\n Switching to use binary mask without potential requirement.')
            allowed_r0_mask = BW;
        end
    end

    % Compute initial points
    if DEBUG
        [r0, ComputeInitialPointsInfo] = computeInitialPoints(POINT_SELECTION_METHOD, allowed_r0_mask, initPointOptions); % Note the points r0 do not consider the mask padding.
    else
        r0 = computeInitialPoints(POINT_SELECTION_METHOD, allowed_r0_mask, initPointOptions); % Note the points r0 do not consider the mask padding.
    end % if

    % If there was only one initial particle, then we do not symulate it,
    % we will just return the centroid of the object
    if size(r0,1) < 2
        [seedPoints,Info] = computeCentroid(BW,DEBUG,2);%'lessThan_2_initial_particles');
        if DEBUG
            Info.ComputeInitialPoints = ComputeInitialPointsInfo;
        end % if
        return;
    end % if

    if DEBUG
        [r_final, Info] = modelParticleDynamics(V, r0, options);

        Info.r0 = r0;
        Info.r_final = r_final;
        Info.V = V;
        Info.ComputeInitialPoints = ComputeInitialPointsInfo;
        Info.message = 0;%'';
    else
        r_final = modelParticleDynamics(V, r0, options);
    end % if

    seedPoints = extractClusterCenters(r_final, options);

catch ME
    % It can be somtimes that there is an error due to a particle flying out of the image region, or something that rarely (hopefully) happens; therefore, we will try to run the function again. If there is an error a second time, then we issue a warning, save the error, and continue on to the next object.
    if errorCount < 1
        [seedPoints, Info] = computeObjectSeedPoints(BW, M, r0set, useCentroid, options, objNumber, errorCount+1);
    else
        seedPoints = [NaN, NaN];
        Info.error = ME;
        if DEBUG
            Info.r0 = [NaN, NaN];
            Info.r_final = [NaN, NaN];
            Info.V = [];
            Info.solverTime = NaN;
            Info.message = 3;%'error';
        end
        fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', objNumber, objNumber)
        fprintf(2,'%s\n', getReport(Info.error,'basic'))
        fprintf('\n')
    end % if
end % try/catch

end % function


function [seedPoint, Info] = computeCentroid(BW,DEBUG,reason)
[i,j] = find(BW);
seedPoint = mean([i,j],1);
Info = [];
if DEBUG
    Info.r0 = [NaN, NaN];
    Info.r_final = [NaN, NaN];
    Info.V = [];
    Info.solverTime = NaN;
    Info.message = reason;
end % if
end % function
