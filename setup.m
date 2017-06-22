function setup

if verLessThan('matlab','9.1')
    error('ceclumpNuclei:setup','The declumpNuclei code requires at least Matlab 2016b because implicit expansion is heavily used.')
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
    error('declumpNuclei:setup','There was an error compiling the required C functions ''interp2mex.c'' and ''nakeinterp1.c''. Make sure that the function ''mex'' is coorectly setup to compile C code.')
end

% Copy histcountsmex into the utilities folder ---------------------------
folder = fullfile(path,'utilities');
d = dir(folder);
names = {d.name};
if ~any(strncmp('DN_histcountsmex',names,16))
    folder = fullfile(matlabroot,'toolbox','matlab','datafun','private');
    d = dir(folder);
    names = {d.name};
    nameIdx = strncmp('histcountsmex',names,13);

    if ~any(nameIdx)
        error('declumpNuclei:setup','File ''histcountsmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the utilities folder and rename it to ''DN_histcountsmex.(extension)''.',folder)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'utilities',['DN_' names{nameIdx}]))
end
fprintf('...added histcountsmex\n')

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
        error('declumpNuclei:setup','File ''pdistmex'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the utilities folder and rename it to ''DN_pdistmex.(extension)''.',folder)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'utilities',['DN_' names{nameIdx}]))
end
fprintf('...added pdistmex\n')

% Add subfolders of current location to path -----------------------------
addpath(genpath(path));
fprintf('...added subfolders to path\n')

fprintf('Setup finished!\n')
    
end
