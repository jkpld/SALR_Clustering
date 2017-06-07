function publishExample(file,reEvaluate)

if nargin < 2
    reEvaluate = true;
end

% Documentation folder
% out_folder = [filesep, 'docs', filesep, 'examples'];
out_folder_img = 'K:\feeling-responsive-gh-pages\images';
out_folder_file = ['K:\feeling-responsive-gh-pages\_posts\example\', datestr(now,29), '-'];

% Get the path of the given file
fileLocation = which(file);

% Geth the path too the given file
pth = string(fileLocation).split(filesep);
pth(end) = [];
pth = char(pth.join(filesep));

% Read in the file text
fid = fopen(fileLocation);
C = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
C = C{1};

% Find the <INSERT FIGURE> tags
fig_locs = find(startsWith(C,'% <INSERT FIGURE'));

% If there are any <INSERT_FIGURE> tags, then we need to replace them with
% code to export the figures and then run the script to create the figures.
if ~isempty(fig_locs) && reEvaluate
    % Copy the file
    file_copy = [pth, filesep, file, '_tempCopy.m'];
    status = copyfile(fileLocation, file_copy);
    if status ~= 1
        error('publishExample:copyFailed','Error publishing file.')
    end

    % Determine how to export the figures.
    %   Use export_fig() if installed (https://github.com/altmany/export_fig)
    %   Otherwise use print()
    if exist('export_fig','file') == 2
        print_fun = @(output) ['export_fig(gcf, ''-png'',''-r200'',''-transparent'', ''', output, ''');'];
    else
        print_fun = @(output) ['print(gcf, ', output, '''-dpng'',''-r200'');'];
    end

    % Replace the lines with the INSER_FIGURE tag with the figure export
    % command.
    C_copy = C;
    fig_name = @(ii) [out_folder_img, filesep, file, sprintf('_%d.png', ii)];
    for ii = 1:length(fig_locs)
        C_copy{fig_locs(ii)} = print_fun(fig_name(ii));
    end

    % Write the new text to the copied file.
    fid = fopen(file_copy,'w');
    cellfun(@(x) fprintf(fid,'%s\n',x),C_copy);
    fclose(fid);

    % Run the file - this should generate the figures we need.
    evalc('run(file_copy)')
    
    % Delete the file copy
    delete(file_copy)
end

% Remove all empty lines at the start of the file.
emptyLine = cellfun(@isempty,C);
if emptyLine(1) == 1
    idx = find(~emptyLine,1,'first');
    C(1:idx-1) = [];
end

% Parse file into code sections
section_heads = startsWith(C,'%%');
figures = startsWith(C,'% <INSERT FIGURE');
thumb_fig = startsWith(C,'% <INSERT FIGURE, THUMB');
comments = startsWith(C,'%') & ~section_heads & ~figures;

line_type = [section_heads, comments, figures, emptyLine];
line_type(:,end+1) = ~any(line_type,2);
line_type = line_type*(1:size(line_type,2))';

header_idx = find(section_heads,1,'first');
C(1:header_idx-1) = [];
line_type(1:header_idx-1) = [];

sections = struct('title',[],'info',[],'content',[],'figure',{});

fig_name = @(ii) [filesep, 'images', filesep, file, sprintf('_%d.png', ii)];

section_counter = 1;
figure_counter = 0;
while 1
    sections(section_counter).title = C(1); % Feed        
    line_type(1) = []; % Eat
    C(1) = []; % Eat

    if isempty(line_type), break; end

    % If next line is comment, then it is part of the sections info.
    if line_type(1) == 2
        end_idx = find(line_type~=2,1,'first');
        if isempty(end_idx)
            end_idx = numel(line_type) + 1;
        end
        sections(section_counter).info = C(1:end_idx-1); % Feed
        line_type(1:end_idx-1) = []; % Eat
        C(1:end_idx-1) = []; % Eat
    end

    if isempty(line_type), break; end

    % Everything until the next section head is content
    end_idx = find(line_type==1,1,'first');
    if isempty(end_idx)
        end_idx = numel(line_type) + 1;
    end

    % Add figures to section
    section_figures = line_type(1:end_idx-1) == 3;
    for ii = sum(section_figures):-1:1
        sections(section_counter).figure{ii} = fig_name(figure_counter + ii);
    end
    figure_counter = figure_counter + sum(section_figures);

    % Remove figure lines from content and save content to section
    sections(section_counter).content = C(~section_figures); % Feed

    line_type(1:end_idx-1) = []; % Eat
    C(1:end_idx-1) = []; % Eat

    if isempty(line_type), break; end
    section_counter = section_counter + 1;
end


if isempty(sections)
    error('publishExample:FoundNoSection', 'Found no sections in the code.')
end


figures = cumsum(figures);
thumb_fig_num = figures(thumb_fig);
if ~isempty(thumb_fig_num)
    thumb_fig_name = [file, sprintf('_%d.png', thumb_fig_num(1))];
else
    thumb_fig_name = '';
end

% Write example
file_out = [out_folder_file, file, '.md'];
fid = fopen(file_out,'w');
try

    write(fid, jekyll_header(sections(1).title, sections(1).info, thumb_fig_name))
    write(fid, code_section(sections(1)))
    
    for i = 2:length(sections)
        write(fid, section_heading(sections(i)))
        write(fid, code_section(sections(i)))
    end

fclose(fid);

catch ME
    fclose(fid);
    rethrow(ME)
end



% sections

% ind = cellfun(@(x) ~isempty(startsWith(x,'% <INSERT FIGURE>')),C);

% if isempty(fileLocation)
%     error('publishExample:file_not_found','The file given, %s, was not found.', file_name)
% end
% 
% pth = string(fileLocation).split(filesep);
% pth(end) = [];
% fileLocation = char(pth.join(filesep))
% 
% 
% 
% file_name
% format = 'html';
% options = struct('format',format,...
%                  'evalCode',true,...
%                  'showCode',false,...
%                  'createThumbnail',false,...
%                  'outputDir', [fileLocation, '\tmp_files\'],...
%                  'figureSnapMethod', 'print',...
%                  'imageFormat', 'png',...
%                  'codeToEvaluate', file_name);
%              
% file = publish(file_name,options);
% file
% % Get the folder the files are located in
%         folder = fileparts(file)
%         
%         % Get the file infomration from the folder
%         file_data = dir(folder)
%         
%         % Get the names of the files in the folder and remove the first two
%         % which are '.' and '..'.
%         names = cell(length(file_data),1);
%         [names{:}] = deal(file_data.name);
%         names(1:2,:) = []
%         
% delete(file)
% 
% 

        

end

function write(fid, C)
    if ~isempty(C)
        cellfun(@(x) fprintf(fid,'%s\n',x),C);
    end
end

function C = section_heading(section)
    C = {};
    if ~isempty(section.title)
        title = strtrim(section.title{1}(3:end));
        C{1} = sprintf('## %s', title);
    end
    
    if ~isempty(section.info)
        tmp = cellfun(@(x) x(2:end), section.info, 'UniformOutput', false);
        C = [C; tmp; ''];
    end
end

function C = code_section(section)
    C = {};
    stuff = [~is_content_empty(section.content), ~isempty(section.figure)];
    
    if ~any(stuff)    
        return;
    end
    
    if stuff(1)
        content = section.content;
        isEmpty = cellfun(@isempty, content);
        start_idx = 1;
        end_idx = length(content);
        
        if isEmpty(1)
            start_idx = find(~isEmpty,1,'first');
        end
        if isEmpty(end)
            end_idx = find(~isEmpty,1,'last');
        end
        content = content(start_idx:end_idx);
        
        content = [' ';'{% highlight matlab %}'; content; '{% endhighlight %}'; ' '];
    end
    
    if stuff(2)
        for ii = length(section.figure):-1:1
            figure_code{ii} = sprintf('<img src="%s">', section.figure{ii});
        end
        figure_code = [''; figure_code; ''];
    end
    
    if all(stuff)
        C = [{...
            '<div class="row">'; ...
            '<div class="medium-7 columns t30" markdown="1">'; ...
            };
            content;
            {...
            '</div>'; ...
            '<div class="medium-5 columns t30">'; ...
            };
            figure_code;
            {...
            '</div>'; ...
            '</div>'}];
    elseif stuff(1)
        C = content;
    else 
        C = figure_code;
    end
end

function C = jekyll_header(title, info, thumb_fig_name)


    title = strtrim(title{1}(3:end));

    C = {...
        '---'; ...
        'layout: code-example'; ...
        sprintf('title: "%s"',title); ...
        'subheadline: "SALR clustering"'; ...
        'categories:'; ...
        '   - example'; ...
        'show_meta: true'};


    if ~isempty(info)
        info = cellfun(@(x) strtrim(x(2:end)), info, 'UniformOutput', false);
        info = string(info).join();
        C = [C; sprintf('teaser: "%s"', info)];
    end
    
    if ~isempty(thumb_fig_name)
        C = [C; 'image:'; sprintf('    thumb: %s', thumb_fig_name)];
    end
    
    C = [C; '---'; ' '];      
end


function out = is_content_empty(C)
    out = false;
    if isempty(C)
        out = true;
    elseif length(C) == 1
        if isempty(C{1})
            out = true;
        end
    end
end
















