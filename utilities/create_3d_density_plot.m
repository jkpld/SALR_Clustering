function Parent = create_3d_density_plot(n, dims, isoLvls, varargin)


[n, dims, isoLvls, markers, color_scale, ticks, ticklabels, cmap, axisPlaneOffset,Parent] = parse_inputs(n,dims,isoLvls,varargin{:});


% Project input data to the three dimensions given -------------------
inds = find(~ismember(1:ndims(n), dims));
for i = inds
    n = sum(n,i);
end

n = squeeze(n);
sz = size(n);
n = permute(n,[2,1,3]);



% Compute the isosurfaces and their bounds ----------------------------
isoBounds = [inf(1,3); -inf(1,3)];
toRemove = false(1,length(isoLvls));

for i = length(isoLvls):-1:1
    iso(i) = isosurface(n,isoLvls(i));
    if isempty(iso(i).vertices)
        toRemove(i) = true;
    else
%         iso(i).vertices = iso(i).vertices(:,[2,1,3]);
        isoBounds(1,:) = min(isoBounds(1,:), min(iso(i).vertices,[],1));
        isoBounds(2,:) = max(isoBounds(2,:), max(iso(i).vertices,[],1));
    end
end

% isoBounds = isoBounds(:,[2,1,3]);

if all(toRemove)
    error('create_d3_density_plot:badIsoLevels','There is no data for any of the isolevels given. Please try changing them.')
end
iso(toRemove) = [];
isoLvls(toRemove) = [];



% Create the figure --------------------------------------------------
ax = axes('Parent',Parent,'Position',[0 0 1 1],'Visible','off');

% Setup the axis and plane bounds
lim_offset = axisPlaneOffset * [-1; 1];
plane_offset = [-2; 2];

plane_bounds = isoBounds + plane_offset;

ax.XLim = isoBounds(:,1) + lim_offset;% * [-1;0];
ax.YLim = isoBounds(:,2) + lim_offset;% * [-1;0];
ax.ZLim = isoBounds(:,3) + lim_offset;% * [-1;0];

% Draw the axis planes
drawAxisPlanes(plane_bounds,ticks,ticklabels,ax)

% Draw a marker to help overlap the isosurfaces with the axis/lines when
% putting the figure together.
line(isoBounds(2,1)+2,isoBounds(2,2)+2,ax.ZLim(1),'marker','.','Color','k','Tag','xy','UserData','isoSurfaceMarker','Visible','off')

% Scale the isolevels for getting good colors (maybe)
isoLvls = isoLvls - min(isoLvls);
isoLvls = isoLvls/max(isoLvls);
cols = interp1(linspace(0,1,size(cmap,1)), cmap, color_scale(isoLvls));
% Plot the isosurfaces and their projects -------------------------------
for i = 1:length(iso)

    patch(iso(i),...
        'FaceColor',cols(i,:),...
        'FaceAlpha',0.3,...
        'EdgeColor','none',...
        'AmbientStrength',0.6,...
        'DiffuseStrength',0.1,...
        'BackFaceLighting','lit',...
        'SpecularStrength',0.5,...
        'SpecularExponent',225,...
        'Tag','isoSurface')
    
    plotProjection(iso(i),cols(i,:),sz,ax);
end



% Plot any markers given ------------------------------------------------
if ~isempty(markers)
    for i = 1:numel(markers)
        if ~isempty(markers(i).dat)
            markers(i).dat = markers(i).dat(:,dims([1,2,3]));
            line(markers(i).dat(:,1),markers(i).dat(:,2),markers(i).dat(:,3),markers(i).options,'Parent',ax)
            if markers(i).project
                projectMarkers(markers(i),ax);
            end
        end
    end
end



% Add lighting ----------------------------------------------------------
for i = 1:2
    for j = 1:2
        for k = 1:2
            light('Position',[ax.XLim(i), ax.YLim(j), ax.ZLim(k)],'Style','local','Parent',ax);
        end
    end
end
lighting gouraud

set(findall(ax,'UserData','contourLine'),'EdgeLighting','none')


% Create a listener to move the planes and the projects when the object is
% rotated ----------------------------------------------------------------
addlistener(ax,'View','PostSet',@viewChange);
addlistener(ax,'XDir','PostSet',@xreverseChange);
addlistener(ax,'YDir','PostSet',@yreverseChange);
ax.UserData.Current_State = [0 0 0 0 0];



% Add axis labels (as the dimension number) ----------------------------
text(mean(ax.XLim), ax.YLim(2), ax.ZLim(1), ['D' num2str(dims(1))],'Tag','xz_xy','UserData',[0,0.15])
text(ax.XLim(2), mean(ax.YLim), ax.ZLim(1), ['D' num2str(dims(2))],'Tag','yz_xy','UserData',[0,0.15])
text(ax.XLim(2), ax.YLim(2), mean(ax.ZLim), ['D' num2str(dims(3))],'Tag','zlbl','Rotation',90,'UserData',[0,0.2])



% Store all of the axis handles for easy use by viewChange() ------------
xy = findall(ax,'Tag','xy');
yz = findall(ax,'Tag','yz');
xz = findall(ax,'Tag','xz');
xz_xy = findall(ax,'Tag','xz_xy');
yz_xy = findall(ax,'Tag','yz_xy');
zlbl = findall(ax,'Tag','zlbl');

change_with_xy = [xy; xz_xy; yz_xy];
change_with_xz = [xz; xz_xy];
change_with_yz = [yz; yz_xy];
change_with_zlbl = zlbl;

ax.UserData.change_with_xy = change_with_xy;
ax.UserData.change_with_xz = change_with_xz;
ax.UserData.change_with_yz = change_with_yz;
ax.UserData.change_with_zlbl = change_with_zlbl;




% Set the default view ---------------------------------------------------
daspect([1 1 1])
view([135,35])
axis vis3d

end

function [n, dims, isoLvls, Markers, ColorScale, Ticks, TickLabels, Colormap,AxisPlaneOffset,Parent] = parse_inputs(n,dims,isoLvls,varargin)

% Validate required inputs -------------------
D = ndims(n);
if length(unique(dims)) ~= 3
    error('There must be at least 3 dimensions given to plot.')
end
if D < 3
    error('Input data should have at least 3 dimensions')
end
validateattributes(isoLvls,{'numeric'},{'nonempty','vector'})


% Setup default colormap ---------------------
if exist('brewermap','file')==2
    cmap = brewermap(9,'GnBu');
    cmap(1,:) = [];
else
    cmap = flipud(parula(9));
end
if all(isoLvls<=1)
    cmap = flipud(cmap);
end

% Setup default ticks --------------------------
sz = size(n);
Ticks = cell(1,3);
TickLabels = cell(1,3);
for i = 1:3
    Ticks{i} = 0:10:sz(i);
    TickLabels{i} = [];
end

% Parse extra inputs ---------------------------
p = inputParser;
p.FunctionName = 'create_3d_density_plot';

addParameter(p,'Markers',[]);
addParameter(p,'ColorScale',@(x) x, @(t) isempty(t) || isa(t,'function_handle'));
addParameter(p,'Ticks',Ticks, @(t) validateattributes(t,{'cell'},{'numel',3}));
addParameter(p,'TickLabels',TickLabels, @(t) validateattributes(t,{'cell'},{'numel',3}));
addParameter(p,'Colormap',cmap, @(t) validateattributes(t,{'double'},{'2d','ncols',3,'real','finite','nonnegative','<=',1}));
addParameter(p,'Parent',[], @(t) isempty(t) || isa(t,'matlab.UI.Figure') || isa(t,'matlab.ui.container.internal.UIContainer') || isa(t,'matlab.ui.container.internal.UIFlowContainer'))
% addParameter(p,'AxisPlaneOffset',10, @(t) validateattributes(t,{'double'},{'numel',1,'real','finite','nonnegative'}))

parse(p,varargin{:})

Markers = p.Results.Markers;
ColorScale = p.Results.ColorScale;
Ticks = p.Results.Ticks;
TickLabels = p.Results.TickLabels;
Colormap = p.Results.Colormap;
AxisPlaneOffset = 10;%p.Results.AxisPlaneOffset;
Parent = p.Results.Parent;

if isempty(Parent)
    Parent = figure('color','w');
    Parent.Units = 'inch';
end
end

function xreverseChange(~,evt)
ax = evt.AffectedObject;

h = ax.UserData.change_with_yz;

for ii = 1:numel(h)
    switch h(ii).Type
        case {'line','patch'}
            currentValue = h(ii).XData(1);
            idx = currentValue ~= ax.XLim;

            h(ii).XData = ax.XLim(idx)*ones(size(h(ii).XData));
        case 'text'
            if ~isempty(h(ii).UserData)
                dr = h(ii).UserData(1);
                dr = dr*[1,-1];
                currentValue = h(ii).Position(1);
                idx = currentValue ~= (ax.XLim + dr);  
                h(ii).Position(1) = ax.XLim(idx) + dr(idx);
            else
                currentValue = h(ii).Position(1);
                idx = currentValue ~= ax.XLim;
                h(ii).Position(1) = ax.XLim(idx);
            end
    end
end

viewChange([],evt)

end

function yreverseChange(~,evt)
ax = evt.AffectedObject;

h = ax.UserData.change_with_xz;

for ii = 1:numel(h)
    switch h(ii).Type
        case {'line','patch'}
            currentValue = h(ii).YData(1);
            idx = currentValue ~= ax.YLim;
            h(ii).YData = ax.YLim(idx)*ones(size(h(ii).YData));
        case 'text'
            if ~isempty(h(ii).UserData)
                dr = h(ii).UserData(1);
                dr = dr*[1,-1];
                currentValue = h(ii).Position(2);
                idx = currentValue ~= (ax.YLim + dr);  
                h(ii).Position(2) = ax.YLim(idx) + dr(idx);
            else
                currentValue = h(ii).Position(2);
                idx = currentValue ~= ax.YLim;
                h(ii).Position(2) = ax.YLim(idx);
            end

    end
end

viewChange([],evt)

end

function viewChange(~,evt)

ax = evt.AffectedObject;

% YZ axis plane
tmp =  mod(abs(floor(ax.View(1)/180)),2) ~= 0;
if tmp ~= ax.UserData.Current_State(1)
    h = ax.UserData.change_with_yz;

    for ii = 1:numel(h)
        switch h(ii).Type
            case {'line','patch'}
                currentValue = h(ii).XData(1);
                idx = currentValue ~= ax.XLim;
                
                h(ii).XData = ax.XLim(idx)*ones(size(h(ii).XData));
            case 'text'
                if ~isempty(h(ii).UserData)
                    dr = h(ii).UserData(1);
                    dr = dr*[1,-1];
                    currentValue = h(ii).Position(1);
                    idx = currentValue ~= (ax.XLim + dr);  
                    h(ii).Position(1) = ax.XLim(idx) + dr(idx);
                else
                    currentValue = h(ii).Position(1);
                    idx = currentValue ~= ax.XLim;
                    h(ii).Position(1) = ax.XLim(idx);
                end
                
        end
    end
end
ax.UserData.Current_State(1) = tmp; % Update value

% XZ axis plane
tmp =  mod(abs(floor((ax.View(1)-90)/180)),2) ~= 0;
if tmp ~= ax.UserData.Current_State(2)
    h = ax.UserData.change_with_xz;

    for ii = 1:numel(h)
        switch h(ii).Type
            case {'line','patch'}
                currentValue = h(ii).YData(1);
                idx = currentValue ~= ax.YLim;
                h(ii).YData = ax.YLim(idx)*ones(size(h(ii).YData));
            case 'text'
                if ~isempty(h(ii).UserData)
                    dr = h(ii).UserData(1);
                    dr = dr*[1,-1];
                    currentValue = h(ii).Position(2);
                    idx = currentValue ~= (ax.YLim + dr);  
                    h(ii).Position(2) = ax.YLim(idx) + dr(idx);
                else
                    currentValue = h(ii).Position(2);
                    idx = currentValue ~= ax.YLim;
                    h(ii).Position(2) = ax.YLim(idx);
                end
                
        end
    end
end
ax.UserData.Current_State(2) = tmp; % Update value


tmp =  [mod(abs(floor((ax.View(1))/90)),2) ~= 0, ...
        mod(abs(floor((ax.View(1))/180)),2) ~= 0];

xlims = ax.XLim;
xf = 1;
if strcmp(ax.XDir,'reverse')
    xlims = flip(xlims);
    xf = -1;
end
    
ylims = ax.YLim;
yf = 1;
if strcmp(ax.YDir,'reverse')
    ylims = flip(ylims);
    yf = -1;
end

h = ax.UserData.change_with_zlbl;
for ii = 1:numel(h)
    switch h(ii).Type
        case 'text'
            dr1 = h(ii).UserData(2);
            dr2 = h(ii).UserData(1);
            dx = dr1*diff(ax.XLim);
            dy = dr1*diff(ax.YLim);

            if isequal(tmp,[0,0])
                h(ii).Position(1:2) = [xlims(1) - xf*dx*cosd(mod((ax.View(1)),90)), ylims(1) + yf*dr2];
            elseif isequal(tmp,[0,1])
                h(ii).Position(1:2) = [xlims(2) + xf*dx*cosd(mod((ax.View(1)),90)), ylims(2) - yf*dr2];
            elseif isequal(tmp,[1,0])
                h(ii).Position(1:2) = [xlims(2) - xf*dr2, ylims(1) - yf*dy*cosd(mod((ax.View(1)),90))];
            elseif isequal(tmp,[1,1])                    
                h(ii).Position(1:2) = [xlims(1) + xf*dr2, ylims(2) + yf*dy*cosd(mod((ax.View(1)),90))];
            end
    end
end


% XY axis plane
tmp = ax.View(2) <= 0;
if tmp ~= ax.UserData.Current_State(3)
    h = ax.UserData.change_with_xy;

    for ii = 1:numel(h)
        switch h(ii).Type
            case {'line','patch'}
                h(ii).ZData = ax.ZLim(tmp+1)*ones(size(h(ii).YData));
            case 'text'
                h(ii).Position(3) = ax.ZLim(tmp+1);
        end
    end
end
ax.UserData.Current_State(3) = tmp; % Update value

h = ax.UserData.change_with_xy;
tmp = ax.View(2) <= 0;
ang = (1-min(1,abs(ax.View(2)/30)))*sign(ax.View(2));
for ii = 1:numel(h)
    if strcmp(h(ii).Type,'text')
        if ~isempty(h(ii).UserData)
            dr = h(ii).UserData(2);
        else
            dr = 0.05;
        end
        h(ii).Position(3) = ax.ZLim(tmp+1) - ang*dr*diff(ax.ZLim);
    end
end

end

function drawAxisPlanes(bnds,ticks,ticklabels,ax)

s = [0 0;
    0 1;
    1 1;
    1 0;
    0 0];

for i = 3:-1:1
    if isempty(ticklabels{i})
        tks{i} = ticks{i};
    else
        tks{i} = ticklabels{i};
    end
    if all(tks{i} == round(tks{i}))
        fmt{i} = '%2.0f';
    else
        fmt{i} = '%2.1f';
    end
    if i~=3
        fmt{i}(2)='0';
    end
end

for dims = 1:3
    
    switch dims
        case 1
            dx = diff(bnds(:,1));
            dz = diff(bnds(:,3));
            z = s(:,1)*dz + min(bnds(:,3));
            x = s(:,2)*dx + min(bnds(:,1));
            y = ax.YLim(1)*ones(size(x));
            
            tag = 'xz';
            
                
            for i = 1:numel(ticks{3})
                if ticks{3}(i) < max(bnds(:,3)) && ticks{3}(i) > min(bnds(:,3))
                    tx = bnds(:,1);
                    ty = ax.YLim(1)*[1;1];
                    tz = ticks{3}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(max(bnds(:,1))+[0;1], ty, tz,'Color',0*[1 1 1],'LineWidth',1)
                    
                    tcklbl = sprintf(fmt{3},tks{3}(i));
                    text(tx(1),ty(1),ticks{3}(i),tcklbl,'Tag','zlbl','UserData',[4,0.1])
                end
            end
            for i = 1:numel(ticks{1})
                if ticks{1}(i) < max(bnds(:,1)) && ticks{1}(i) > min(bnds(:,1))
                    tz = bnds(:,3);
                    ty = ax.YLim(1)*[1;1];
                    tx = ticks{1}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(tx, ty, max(bnds(:,3))+[0;1],'Color',0*[1 1 1],'LineWidth',1)
                end
            end
            
        case 2
            
            dy = diff(bnds(:,2));
            dz = diff(bnds(:,3));
            z = s(:,1)*dz + min(bnds(:,3));
            y = s(:,2)*dy + min(bnds(:,2));
            x = ax.XLim(1)*ones(size(y));            
            
            tag = 'yz';
            for i = 1:numel(ticks{3})
                
                if ticks{3}(i) < max(bnds(:,3)) && ticks{3}(i) > min(bnds(:,3))
                    ty = bnds(:,2);
                    tx = ax.XLim(1)*[1;1];
                    tz = ticks{3}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(tx, max(bnds(:,2))+[0;1], tz,'Color',0*[1 1 1],'LineWidth',1)
                end
            end
            for i = 1:numel(ticks{2})
                if ticks{2}(i) < max(bnds(:,2)) && ticks{2}(i) > min(bnds(:,2))
                    tz = bnds(:,3);
                    tx = ax.XLim(1)*[1;1];
                    ty = ticks{2}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(tx, max(bnds(:,2))+[0;1], tz,'Color',0*[1 1 1],'LineWidth',1)
                end
            end
            
        case 3
            
            dx = diff(bnds(:,1));
            x = s(:,1)*dx + min(bnds(:,1));
            
            dy = diff(bnds(:,2));
            y = s(:,2)*dy + min(bnds(:,2));
            
            z = ax.ZLim(1)*ones(size(y));
            tag = 'xy';
            for i = 1:numel(ticks{1})
                
                if ticks{1}(i) < max(bnds(:,1)) && ticks{1}(i) > min(bnds(:,1))
                    ty = bnds(:,2);
                    tz = ax.ZLim(1)*[1;1];
                    tx = ticks{1}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(tx, max(bnds(:,2))+[0;1], tz,'Color',0*[1 1 1],'LineWidth',1)
                    tcklbl = sprintf(fmt{1},tks{1}(i));
                    text(tx(1),ax.YLim(2)-5,tz(1),tcklbl,'Tag','xz_xy','UserData',[5,0.05],'HorizontalAlignment','center')
                end
            end
            for i = 1:numel(ticks{2})   
                if ticks{2}(i) < max(bnds(:,2)) && ticks{2}(i) > min(bnds(:,2))
                    tx = bnds(:,1);
                    tz = ax.ZLim(1)*[1;1];
                    ty = ticks{2}(i)*[1;1];
                    
                    line(tx, ty, tz,'Color',0.8*[1 1 1],'LineWidth',0.5,'Tag',tag)
                    % line(max(bnds(:,1))+[0;1], ty, tz,'Color',0*[1 1 1],'LineWidth',1)
                    tcklbl = sprintf(fmt{2},tks{2}(i));
                    text(ax.XLim(2)-5,ty(1),tz(1),tcklbl,'Tag','yz_xy','UserData',[5,0.05],'HorizontalAlignment','center')
                end
            end
    end
    
    line(x,y,z,'Color','k','LineWidth',1.5,'Parent',ax,'Tag',tag)
    
    %     patch('Vertices',[x(1:end-1),y(1:end-1),z(1:end-1)],'Faces',[1,2,3,4],'FaceColor',0.99*[1 1 1],'EdgeColor','k','LineWidth',1.5,'Parent',ax,'FaceLighting','none','EdgeLighting','none')
end

end

function plotProjection(fv,col,sz,ax)

v = round(fv.vertices);

for ind = 1:3

    inds = find(~ismember(1:3, ind));

    tmp = zeros(sz(inds));

    tmp(v(:,inds(1)) + (v(:,inds(2))-1)*sz(inds(1))) = 1;
    
    cntrs = bwboundaries(tmp);
    
    kappaSmoothingSigma=0.7;
    filtSize = round(7*kappaSmoothingSigma);
    
    x = -floor(filtSize/2):ceil(filtSize/2);
    G = exp(-(x).^2/(2*kappaSmoothingSigma^2))/(kappaSmoothingSigma*sqrt(2*pi));
    G = G(:);
    
    for c = 1:numel(cntrs)
        cntr = cntrs{c};

        % Get boundary curvature
        cntr = imfilter(cntr,G,'circular','conv','same');
        cntr = [cntr;cntr(1,:)]';


        switch ind
            case 2
                z = cntr(2,:);
                x = cntr(1,:);
                y = ax.YLim(1)*ones(size(x));
                tag = 'xz';

            case 1
                z = cntr(2,:);
                y = cntr(1,:);
                x = ax.XLim(1)*ones(size(y));
                tag = 'yz';
            case 3
                x = cntr(1,:);
                y = cntr(2,:);
                z = ax.ZLim(1)*ones(size(x));
                tag = 'xy';
        end
        patch('XData',x,'YData',y,'ZData',z,'FaceColor','none','EdgeColor',col,'Tag',tag,'Parent',ax,'UserData','contourLine','LineWidth',2,'EdgeLighting','none','FaceLighting','none')
%         line(x,y,z,'Color',col,'LineWidth',2,'Tag',tag,'Parent',ax,'UserData','contourLine')
    end
end
end

function projectMarkers(markers,ax)
dat = markers.dat;
options = markers.options;

for i = 1:3
    switch i
        case 1
            y = dat(:,2);
            z = dat(:,3);
            x = ax.XLim(1)*ones(size(y));
            tag = 'yz';
        case 2
            x = dat(:,1);
            z = dat(:,3);
            y = ax.YLim(1)*ones(size(x));
            tag = 'xz';
        case 3
            x = dat(:,1);
            y = dat(:,2);
            z = ax.ZLim(1)*ones(size(x));
            tag = 'xy';
    end
    line(x,y,z,options,'Parent',ax,'Tag',tag)
end
end

function plotCntrLines(n,dims,lvls,cmap,ax)

numC = size(cmap,1);

inds = find(~ismember(1:ndims(n), dims));
for i = inds
    n = sum(n,i);
end
n = squeeze(n);

C = contourc(double(n), lvls);%

cntrs = {};
cols = {};

counter = 1;
while true
    lvl = C(1,1);
    N = C(2,1);
    
    x = C(1,2:N+1);
    y = C(2,2:N+1);
    C(:,1:N+1) = [];
    
    cols{counter} = interp1(linspace(0,1,numC),cmap,(lvl/max(lvls))^(1/2));
    cntrs{counter} = [x;y];
    
    counter = counter + 1;
    
    if isempty(C)
        break;
    end
end

cntrs = flip(cntrs);
cols = flip(cols);


for i = 1:numel(cntrs)
    switch inds
        case 1
            z = cntrs{i}(1,:);
            x = cntrs{i}(2,:);
            y = ax.YLim(1)*ones(size(x));
            tag = 'xz';
            
        case 2
            z = cntrs{i}(1,:);
            y = cntrs{i}(2,:);
            x = ax.XLim(1)*ones(size(y));
            tag = 'yz';
        case 3
            x = cntrs{i}(1,:);
            y = cntrs{i}(2,:);
            z = ax.ZLim(1)*ones(size(x));
            tag = 'xy';
    end
    line(x,y,z,'Color',cols{i},'LineWidth',2,'Tag',tag)
end

end