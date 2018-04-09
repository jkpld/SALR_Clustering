function setup
% SETUP Add required files to the path and try to compile .c functions for
% increased speed.
%
% setup()

% James Kapaldo

% Check dependancies -----------------------------------------------------
if verLessThan('matlab','9.1')
    error('SALR_Clustering:setup','The SALR_Clustering code requires at least Matlab 2016b because implicit expansion is heavily used.')
end
verImg = ver('images');
if isempty(verImg)
    error('SALR_Clustering:setup','The Image Processing Toolbox is required but not found.')
end
verStats = ver('stats');
if isempty(verStats)
    error('SALR_Clustering:setup','The Statistics and Machine Learning Toolbox is required but not found.')
end

fprintf('Setting up...\n')

% Get the path to the current folder -------------------------------------
fileLocation = mfilename('fullpath');
path = fileparts(fileLocation);

fprintf('...Found path\n')

% Compile the needed C files into .mex files -----------------------------
have_compiler = true;
try 
    evalc('mex(''-setup'',''c'')');
catch % ME
    have_compiler = false;
    warning('SALR_Clustering:setup','There is no supported C compiler installed.\nThe code will still run; however, it could be slower for 2D data without the compiled mex functions.')
%     rethrow(ME)
end

if have_compiler
    try
%         if ~exist(fullfile(path,'utilities',['interp2mex.' mexext]),'file')
            mex(fullfile(path,'utilities','interp2mex.c'),'-outdir',fullfile(path,'utilities'),'-silent')
%         end
        fprintf('...Compiled interp2mex.c\n')
%         if ~exist(fullfile(path,'utilities',['nakeinterp1.' mexext]),'file')
            mex(fullfile(path,'utilities','nakeinterp1.c'),'-outdir',fullfile(path,'utilities'),'-silent')
%         end
        fprintf('...Compiled nakeinterp1.c\n')
    catch % ME
    %     rethrow(ME)
        warning('SALR_Clustering:setup','There was an error compiling the required C functions ''interp2mex.c'' and ''nakeinterp1.c''. Make sure that the function ''mex'' is coorectly setup to compile C code.\nThe code will still run; however, it could be slower for 2D data without the compiled mex functions.')
    end
end

% Copy pdistmex into the utilities folder --------------------------------
folder = fullfile(path,'utilities');
d = dir(folder);
names = {d.name};
% if ~any(strncmp('salr_pdistmex',names,11))
    folder = fullfile(matlabroot,'toolbox','stats','stats','private');
    d = dir(folder);
    names = {d.name};
    nameIdx = strncmp('pdistmex',names,8);

    if ~any(nameIdx)
        error('SALR_Clustering:setup','File ''pdistmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the utilities folder and rename it to ''salr_pdistmex.%s''.',folder,mexext)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'utilities',['salr_' names{nameIdx}]))
% end
fprintf('...Added pdistmex\n')

% Add subfolders of current location to path -----------------------------
% (but do not include any .git repositories
pths = split(string(genpath(path)),';');
toIgnore = {'.git','docs'};
pths = pths(~pths.contains(toIgnore)).join(';');
addpath(pths.char());
fprintf('...Added subfolders to path\n')

fprintf('Setup finished!\n')
    
end
