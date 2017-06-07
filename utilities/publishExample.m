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
    evalc('run(file_copy)');

    % Clean up
    delete(file_copy)
    close all
end

% Parse file into code sections -----------------------------------------
% Get the empty lines
emptyLine = cellfun(@isempty,C);
% Section heads
section_heads = startsWith(C,'%%');
% figure locations
figures = startsWith(C,'% <INSERT FIGURE');
thumb_fig = startsWith(C,'% <INSERT FIGURE, THUMB'); % figure to use for thumbnail image

% Captions
caption_starts = find(startsWith(C,'% <CAPTION>'));
caption_ends = find(startsWith(C,'% </CAPTION>'));

if numel(caption_starts) ~= numel(caption_ends)
    error('publishExample:BadCaption','The number of caption start tags <CAPTION>, %d, is not the same as the number of caption end tags </CAPTION>, %d.', numel(caption_starts), numel(caption_ends))
end
if caption_ends(1) < caption_starts(1)
    error('publishExample:BadCaption','A caption end </CAPTION> appears before a caption start <CAPTION>.')
end
captions = false(numel(C),1);
for ii = 1:numel(caption_starts)
    captions(caption_starts(ii):caption_ends(ii)) = true;
end

% Comments
comments = startsWith(C,'%') & ~section_heads & ~figures & ~captions;

line_type = [section_heads, comments, figures, captions, emptyLine];
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
    section_figures_idx = find(section_figures);
    for ii = numel(section_figures_idx):-1:1
        % Does the figure have a caption?
        if line_type(section_figures_idx+1) == 4
            cap_end_idx = section_figures_idx + find(line_type(section_figures_idx+1:end)~=4,1,'first') - 1;
            if isempty(cap_end_idx)
                cap_end_idx = numel(line_type);
            end
            caption = C(section_figures_idx+2:cap_end_idx-1);
            caption = cellfun(@(x) strtrim(x(2:end)), caption, 'UniformOutput', 0);
            caption = char(string(caption).join(' '));
        else
            caption = '';
        end
        sections(section_counter).figure{ii} = {fig_name(figure_counter + ii), caption};
    end

    figure_counter = figure_counter + numel(section_figures_idx);

    % Remove figure/caption lines from content and save content to section
    section_captions = line_type(1:end_idx-1) == 4;
    sections(section_counter).content = C(~section_figures & ~section_captions); % Feed

    line_type(1:end_idx-1) = []; % Eat
    C(1:end_idx-1) = []; % Eat

    if isempty(line_type), break; end
    section_counter = section_counter + 1;
end
% error('some error')
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
        if ~isempty(title)
            C{1} = sprintf('## %s', title);
        end
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
        toIndent = startsWith(content,'''');
        content(toIndent) = cellfun(@(x) ['    ' x], content(toIndent),'UniformOutput',false);

        content = [' ';'{% highlight matlab %}'; content; '{% endhighlight %}'; ' '];
    end

    if stuff(2)
        figure_code = {};
        for ii = length(section.figure):-1:1
            figure_code = [figure_code; {sprintf('<img src="%s">', section.figure{ii}{1})}];

            if ~isempty(section.figure{ii}{2})
                cap_code = {'<figcaption class="text-right">'; ...
                    section.figure{ii}{2}; ...
                    '</figcaption>'};
                figure_code = [figure_code; cap_code];
            end
        end
        figure_code = [' '; figure_code; ' '];
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
