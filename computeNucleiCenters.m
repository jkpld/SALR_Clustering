function [seedPoints,Info] = computeNucleiCenters(I,BW,options)
% COMPUTENUCLEICENTERS Declump the nuclei in an image
%
% [BW,cuts,Info] = computeNucleiCenters(I,BW,options)
%
% I - input image
% BW - object mask for image
% options - declumpOptions class object
%
% BW - object mask after partitioning clumps
% cuts - Nx4 array. cuts(i,1:2) and cuts(i,3:4) give the two vertices of
%        the i'th cut
% Info - cell array of structures giving information

% James Kapaldo

% if options.Use_GPU
%     dev = gpuDevice();
% end

% Get number of rows in image
numImRows = size(I,1);

% Smooth the image
Is = imfilter(I,fspecial('gaussian',7,1));

% Remove small holes from mask
BW = ~bwareaopen(~BW,options.Minimum_Hole_Size,4);
CC = bwconncomp(BW); % Get connected commponents
pixelList = CC.PixelIdxList;

% Create interior confining potential -----------------------------------
V_interior = create_base_interior_confining_potential(BW); % Base potential
H = homogenize_image(Is,BW,CC,options); % Get homogenized image
V_interior = V_interior ./ H;

% Create a sliced cell array of V_interior for each object
V_interior = cellfun(@(x) V_interior(x), pixelList,'UniformOutput',false);

% Compute boundary information
[B,~,K,r0set] = computeBoundaryInformation(BW,options);


% If there is an object of interest, remove all the others.
if ~isempty(options.Object_Of_Interest)
    B = B(options.Object_Of_Interest);
    K = K(options.Object_Of_Interest);
    r0set = r0set(options.Object_Of_Interest);
    pixelList = pixelList(options.Object_Of_Interest);
    options.Use_Parallel = false;
end

% Get object area and determine if the object is convex (no positive
% curvature)
area = cellfun(@numel,pixelList);
isConvex = cellfun(@(x) ~any(x>0),K);
useCentroid = isConvex || (area < pi*options.Wigner_Seitz_Radius.^2);

% Get object offset
objOffset = cellfun(@(x) min(B,[],1,'omitnan') - 1, B,'UniformOutput',false);

% Offset the r0set points to coorespond to object origin and not image
% origin.
r0set = cellfun(@(x,y) x - y, r0set, objOffset,'UniformOutput',false);

% Initialize sliced variables
Info = cell(numel(B),1);
seedPoints = cell(numel(B),1);

% Convert the options to a normal structure so that we dont have to copy it
% and re-intialize each copy.
if options.Use_Parallel
    warning('off','MATLAB:structOnObject')
    options = struct(options);
    warning('on','MATLAB:structOnObject')
end

% Declump objects

if options.Use_Parallel

    parfor obj = 1:numel(B)

        % Create mask and interior potential images for the object
        [objBW,V] = createObjectImages(pixelList{obj}, numImRows, true(numel(pixList),1), V_interior{obj});
            
        if useCentroid(obj)
            [i,j] = find(objBW);
            seedPoints{obj} = [mean([i,j]) + objOffset{obj}, obj];
            Info{obj}.r0 = [];
            Info{obj}.r_end = [];
            Info{obj}.solverTime = NaN;
            Info{obj}.convexOrTooSmall = true;
            Info{obj}.error = [];
        else    
            [seedPoints{obj},Info{obj}] = copmuteObjectSeedPoints(objBW, V, r0set{obj}, options);

            if ~isempty(Info{obj}.error)
                % If there was an error, try again as most error occure
                % because the centers were in just the wrong position and
                % they will not be again.
                [seedPoints{obj},Info{obj}] = copmuteObjectSeedPoints(objBW, V, r0set{obj}, options);

                if ~isempty(Info{obj}.error)
                    fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', obj, obj)
                    fprintf(2,'%s\n', getReport(Info{obj}.error,'basic'))
                    fprintf('\n')
                end % if error
            end % if error
            
            Info{obj}.convexOrTooSmall = false;
            
            % Shift r0, r_end and seed_points back to image coordinates :
            % could also do this with some cellfun()'s after the loop too
            if ~isempty(Info{obj}.r0)
                Info{obj}.r0 = Info{obj}.r0 + objOffset{obj};
                Info{obj}.r_end = Info{obj}.r_end + objOffset{obj};
            end
            if ~isempty(seedPoints{obj})
                seedPoints{obj} = [seedPoints{obj} + objOffset{obj}, obj*size(seedPoints{obj},1)];
            end
        end % if useCentroid
    end % parfor
else

    generateDisplayAt = unique(round(linspace(1,numel(B),7)));
    processTimes = zeros(1,numel(B));

    for obj = 1:numel(B)
        procTime = tic;

        % Create mask and interior potential images for the object
        [objBW,V] = createObjectImages(pixelList{obj}, numImRows, true(numel(pixList),1), V_interior{obj});
            
        if useCentroid(obj)
            [i,j] = find(objBW);
            seedPoints{obj} = mean([i,j]);
            Info{obj}.r0 = [];
            Info{obj}.r_end = [];
            Info{obj}.solverTime = NaN;
            Info{obj}.convexOrTooSmall = true;
            Info{obj}.error = [];
        else    
            [seedPoints{obj},Info{obj}] = copmuteObjectSeedPoints(objBW, V, r0set{obj}, options);

            if ~isempty(Info{obj}.error)
                % If there was an error, try again as most error occure
                % because the centers were in just the wrong position and
                % they will not be again.
                [seedPoints{obj},Info{obj}] = copmuteObjectSeedPoints(objBW, V, r0set{obj}, options);

                if ~isempty(Info{obj}.error)
                    fprintf('\nWarning! There was an error in object %d. Full error report stored in Info{%d}.error\n', obj, obj)
                    fprintf(2,'%s\n', getReport(Info{obj}.error,'basic'))
                    fprintf('\n')
                end % if error
            end % if error
            
            Info{obj}.convexOrTooSmall = false;
            
            % Shift r0, r_end and seed_points back to image coordinates :
            % could also do this with some cellfun()'s after the loop too
            if ~isempty(Info{obj}.r0)
                Info{obj}.r0 = Info{obj}.r0 + objOffset{obj};
                Info{obj}.r_end = Info{obj}.r_end + objOffset{obj};
            end
            if ~isempty(seedPoints{obj})
                seedPoints{obj} = [seedPoints{obj} + objOffset{obj}, obj*size(seedPoints{obj},1)];
            end
        end % if useCentroid

        processTimes(obj) = toc(procTime);
        if any(obj == generateDisplayAt)
            fprintf('%s >> %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),obj,numel(B),sum(processTimes(1:obj))/60,mean(processTimes(1:obj))*numel(B)/60)
        end % if
    end % for
end % if

if options.Use_GPU
    gpuDevice([]);
end

% Produce output
seedPoints = cat(1,seedPoints{:});

if all(cellfun(@isempty, Info))
    Info = [];
end

end



function [seedPoints, Info] = copmuteObjectSeedPoints(BW,V,r0set,options)

try

    % Initialize values
    seedPoints = [];

    if options.Debug
        Info = struct('r0',[],'r_end',[],'solverTime',NaN,'error',[]);
    else
        Info.error = [];
    end

    % TODO: Add both of these parameters fo the declumpOptions (and rename
    % declumpOptions to seedPointOptions)
    PAD_SIZE = options.Potential_Pad_Size;
    R0_MAXV = options.Maximum_r0_Potential;
    
    % Create full confining potential ---------------------------------------
    V = add_exterior_confining_potential(BW,V,options); % Remember that this output V has been padded

    % Get initial particle locations ----------------------------------------
    initPointOptions.rs = options.Wigner_Seitz_Radius;
    initPointOptions.r0set = r0set;

    % Don't let the particles start with too high a potential energy.
    allowed_r0_mask = BW & ( V(PAD_SIZE+1:end-PAD_SIZE,PAD_SIZE+1:end-PAD_SIZE) < R0_MAXV );


    if options.Debug
        [r0, ComputeInitialPointsInfo] = computeInitialPoints(options.Point_Selection_Method, allowed_r0_mask, initPointOptions); % Note the points r0 do not consider the mask padding.
        [seedPoints, Info] = computeSeedPoints(V, r0, options);

        Info.r0 = r0;
        Info.ComputeInitialPoints = ComputeInitialPointsInfo;
        Info.error = [];
    else
        r0 = computeInitialPoints(options.Point_Selection_Method, allowed_r0_mask, initPointOptions); % Note the points r0 do not consider the mask padding.
        seedPoints = computeSeedPoints(V, r0, options);
    end



catch ME
    Info.error = ME;
end

end
