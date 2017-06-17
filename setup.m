function setup
% SETUP Add required files to the path and try to compile .c functions for
% increased speed.
%
% setup()

% James Kapaldo

if verLessThan('matlab','9.1')
    error('seed_point_detection:setup','The seed_point_detection code requires at least Matlab 2016b because implicit expansion is heavily used.')
end

fprintf('\nStarting setup...\n')

% Get the path to the current folder -------------------------------------
fileLocation = mfilename('fullpath');
path = fileparts(fileLocation);

fprintf('...found path\n')

% Compile the needed C files into .mex files -----------------------------
try 
    mex('-setup','c')
catch ME
    fprintf(2,'Error! There is no supported C compiler installed. Please install a supported C compiler to continue.\n\n')
    rethrow(ME)
end

try
    if ~exist(fullfile(path,'utilities',['interp2mex.' mexext]),'file')
        mex(fullfile(path,'utilities','interp2mex.c'),'-outdir',fullfile(path,'utilities'),'-silent')
    end
    fprintf('...compiled interp2mex.c\n')
    if ~exist(fullfile(path,'utilities',['nakeinterp1.' mexext]),'file')
        mex(fullfile(path,'utilities','nakeinterp1.c'),'-outdir',fullfile(path,'utilities'),'-silent')
    end
    fprintf('...compiled nakeinterp1.c\n')
catch ME
%     rethrow(ME)
    error('seed_point_detection:setup','There was an error compiling the required C functions ''interp2mex.c'' and ''nakeinterp1.c''. Make sure that the function ''mex'' is coorectly setup to compile C code.')
end

% Copy pdistmex into the utilities folder --------------------------------
folder = fullfile(path,'utilities');
d = dir(folder);
names = {d.name};
if ~any(strncmp('DN_pdistmex',names,11))
    folder = fullfile(matlabroot,'toolbox','stats','stats','private');
    d = dir(folder);
    names = {d.name};
    nameIdx = strncmp('pdistmex',names,8);

    if ~any(nameIdx)
        error('seed_point_detection:setup','File ''pdistmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the utilities folder and rename it to ''DN_pdistmex.(extension)''.',folder)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'utilities',['DN_' names{nameIdx}]))
end
fprintf('...added pdistmex\n')

% Add subfolders of current location to path -----------------------------
% (but do not include any .git repositories
pths = split(string(genpath(path)),';');
toIgnore = {'.git','docs'};
pths = pths(~pths.contains(toIgnore)).join(';');
addpath(pths.char());
fprintf('...added subfolders to path\n')

fprintf('Setup finished!\n')
    
end
