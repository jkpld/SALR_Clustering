function [seedPoints, Info] = processObjects(pixels, M, r0set, useCentroid, nRows, options)

% See also COMPUTEOBJECTSEEDPOINTS CREATEOBJECTIMAGES MODELPARTICLEDYNAMICS

    % Number of input objects
    N = numel(pixels);

    % If there is only one object, then turn off parallel computing.
    if N == 1
        options.Use_Parallel = false;
    end


    % If we are computing in parallel, then first convert the options class
    % element to a structure to prevent reinitiallization on transfer to
    % each worker.
    if options.Use_Parallel
        warning('off','MATLAB:structOnObject')
        options = struct(options);
        warning('on','MATLAB:structOnObject')
    end

    % Initialize sliced variables
    Info = cell(N,1);
    seedPoints = cell(N,1);

    if options.Use_Parallel
        parfor obj = 1:N
            
            % Create mask and interior potential images for the object
            [objBW,objM] = createObjectImages(pixels{obj}, nRows, true(numel(pixels{obj}),1), M{obj});
            [seedPoints{obj}, Info{obj}] = computeObjectSeedPoints(objBW, objM, r0set{obj}, useCentroid(obj), options, obj)

            % Add object number as third column. (This is mostly just helpful when comparing against truth data, as the truth data is labeled by each object.)
            seedPoints{obj} = [seedPoints{obj}, obj*ones(size(seedPoints{obj},1),1)];

            % ===========================================================
            % This could be a good place to put object segmentation code.
            % ===========================================================
        end
    else
        % Display the progress of the calculation (we can do this since we are not computing the objects in parallel)
        generateDisplayAt = unique(round(linspace(1,N,7)));
        processTimes = zeros(1,N);
        fprintf('Starting seed point calculation\n')

        for obj = 1:N
            procTime = tic;

            % Create mask and interior potential images for the object
            [objBW,objM] = createObjectImages(pixels{obj}, nRows, true(numel(pixels{obj}),1), M{obj});
            [seedPoints{obj}, Info{obj}] = computeObjectSeedPoints(objBW, objM, r0set{obj}, useCentroid(obj), options, obj);

            % Add object number as third column. (This is mostly just helpful when comparing against truth data, as the truth data is labeled by each object.)
            seedPoints{obj} = [seedPoints{obj}, obj*ones(size(seedPoints{obj},1),1)];

            % ===========================================================
            % This could be a good place to put object segmentation code.
            % ===========================================================

            processTimes(obj) = toc(procTime);
            if any(obj == generateDisplayAt)
                fprintf('%s >> %d/%d (%0.2f/%0.2f)...\n',datestr(now,31),obj,N,sum(processTimes(1:obj))/60,mean(processTimes(1:obj))*N/60)
            end % if
        end % for
        fprintf('Finished!\n')
    end % if
end % function
