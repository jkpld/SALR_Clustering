function publishDocumentation
% PUBLISHDOCUMENTATION  Function to create the documentation data used in
% creating the website.
%
% publishDocumentation()
%
% This function is used to create the code.yml file in docs/_data

% James Kapaldo

fileLocation = mfilename('fullpath');
path = fileparts(fileparts(fileLocation));



toInclude = {'coreFunctions','confiningPotential'};
output_pth = fullfile(path, 'docs','_data','code.yml');

fid = fopen(output_pth,'w');
% if fid > -1

try
    writeLine(fid, 'title', 'Getting Started', 0);
    writeLine(fid, 'url', '/quick-start/', 0);
    writeLine(fid, 'title', 'Code', 0);
    writeLine(fid, 'url', '/code/', 0);
    writeLine(fid, 'dropdown', '', 0);

    writeLine(fid, 'title', 'Interface functions', 1);
    writeLine(fid, 'url', '/interface-functions/', 1);
    addDirFunctions(fid, path, 1)

    writeLine(fid, 'title', 'Core functions', 1);
    writeLine(fid, 'url', '/core-functions/', 1);
    addDirFunctions(fid, fullfile(path, 'coreFunctions'), 1)

    writeLine(fid, 'title', 'Setup functions', 1);
    writeLine(fid, 'url', '/setup-functions/', 1);
    addDirFunctions(fid, fullfile(path, 'setupFunctions'), 1)

    writeLine(fid, 'title', 'Utilities', 1);
    writeLine(fid, 'url', '/utility-functions/', 1);
    addDirFunctions(fid, fullfile(path, 'utilities'), 1)

    writeLine(fid, 'title', 'Examples', 0);
    writeLine(fid, 'url', '/examples/', 0);
    
    writeLine(fid, 'title', 'Data documentation', 0);
    writeLine(fid, 'url', '/data-documentation/', 0);


catch ME
    fclose(fid);
    rethrow(ME);
end

fclose(fid);
% else
%     error('publishDocumentation:cantOpen','Failed to open file for writing.');
% end


end

function addDirFunctions(fid, folderName, offset)
D = dir(fullfile(folderName, '*.m'));
names = {D.name}';

if ~isempty(names)
    writeLine(fid, 'functions', '', offset);
    for jj = 1:numel(names)
        file = fullfile(folderName, names{jj});
        str = shortDescription(file);
        C = functionDeclarations(file);
        writeLine(fid, 'title', names{jj}(1:end-2), offset+1);
        writeLine(fid, 'teaser', str, offset+1);
        if ~isempty(C)
            writeLine(fid, 'definitions', '', offset+1);
            str = strjoin(C,'\\n');
            writeLine(fid, 'def', str, offset+2);
        end
    end
end
end

function str = shortDescription(file)

[~,name] = fileparts(file);
str = help(file);
str = strsplit(str,'\n')';
str = cellfun(@(x) strtrim(x),str,'UniformOutput',0);
emptyLines = cellfun(@isempty, str);
str = str(1:find(emptyLines,1,'first')-1);
str = strjoin(str);
str = strrep(str,upper(name),'');
str = strtrim(str);

end

function C = functionDeclarations(file)

%     [~,name] = fileparts(file);
str = help(file);
str = strsplit(str,'\n')';
str = cellfun(@(x) strtrim(x),str,'UniformOutput',0);

% Regular expression for matching lines that call a function.
% Ex: It will match the following three lines
%
%    [a,b] = foo(c)
% c = foo2(a, b);
%  [d] = foo3(a, c, ...)
%
% but it will not match the below lines
%
% as  [a,b] = foo(c)
% = foo2(c)
% [] = foo3()
%    [a,b] = foo(c) lkj

% pattern1 = '^[ ]*(\[\w+.*?\]|\w+)[ ]*=[ ]*\w+(\([\w\.]+?.*?\))\;*[ ]*(|\%.*)$';
pattern1 = '^[ ]*(|((\[\w+.*?\]|\w+)[ ]*=[ ]*))\w+(\([^\)]*?\))\;*[ ]*(|\%.*)$';
% can start with any amount of spaces: ^[ ]*
% has either square brackets enclosing at least one word or a single
%    word: (\[\w+.*?\]|\w+)
% an equal sign surrounded by any number of spaces: [ ]*=[ ]*
% one word: \w+
% parenthesis surrounding at least one word or dots: (\([\w\.]+.*?\))
% an optional semicolon: \;*
% any number of white spaces till end of line: [ ]*$


% Also, match lines that look like a variable description, like
%  y = [ a, b, c ]

% pattern2 = '^[ ]*\w+[ ]*=[ ]*(\[[ ]*\w+.*?\])''*\;*[ ]*(|\%.*)$';
% pattern = sprintf('(%s|%s)', pattern1, pattern2);

% Ingore pattern2 for now ...
matches = regexp(str, pattern1, 'match');
C = str(~cellfun(@isempty,matches));
C = cellfun(@(x) strtrim(x),C,'UniformOutput',0);

end

function writeLine(fid,type,str,level)

line = repmat(' ',1,2*(level+1));

if any(strcmp(type,{'title','def'}))
    line(end-1) = '-';
end
if ~isempty(str)
    str = ['"', str, '"'];
end
fprintf(fid, '%s%s: %s\n', line, type, str);

end


%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
