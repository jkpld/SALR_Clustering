function goDark(fig)


try
    if strcmp(fig.Type,'figure')
        fig.Color = [0 0 0];
    end
catch
end

childs = [findall(fig,'type','text');...
    findall(fig,'Style','text');...
    findall(fig,'type','textbox');...
    findall(fig,'type','axes');...
    findall(fig,'type','colorbar');...
    findall(fig,'type','uiflowcontainer');...
    findall(fig,'type','uicontainer'); ...
    findall(fig,'type','uigridcontainer'); ...
    findall(fig,'type','uicontrol');...
    findall(fig,'type','legend')];

for i = 1:length(childs)
    
    h = childs(i);
    
    if strcmp(h.Tag, 'ignore')
        continue;
    end
    
    try
        switch h.Style
            case 'text'
                h.ForegroundColor = [0.6 0.6 0.6];
                h.BackgroundColor = [0 0 0];
            case 'textbox'
                h.Color = 0.6*[1 1 1];
            case {'radiobutton','checkbox'}
                h.BackgroundColor = 0.1*[1 1 1];
        end
    catch 
    end
    
    try
        switch h.Type
            case 'text'
                h.Color = [0.6 0.6 0.6];
            case {'TextBox', 'textboxshape'}
                h.Color = [0.6 0.6 0.6];
            case 'axes'
                try
                    h.Backdrop.FaceColor = [0.2 0.2 0.2];
                    h.GridColor = [0.5 0.5 0.5];
                    h.YColor = [0.6 0.6 0.6];
                    h.XColor = [0.6 0.6 0.6];
                    h.ZColor = [0.6 0.6 0.6];
                    h.Title.Color = [0.9 0.9 0.9];
                    h.Box = 'off';
                    try
                        h.XRuler.TickDir = 'out';
                        h.YRuler.TickDir = 'out';
                    catch
                        h.TickDir = 'out';
                    end

                    h.XGrid = 'on';
                    h.YGrid = 'on';
                    h.ZGrid = 'on';
                catch
                end
            case 'colorbar'
                
                h.Ruler.Color = [0.6 0.6 0.6];
                h.BoxHandle.Visible = 'off';
                h.TickLength = 0.01;
                try
                    h.TickDirection = 'out';
                catch
                    h.Ruler.TickDir = 'out';
                end
                h.Ruler.Axle.Visible = 'off';
                
            case {'uiflowcontainer','uicontainer','uigridcontainer'}
                h.BackgroundColor = 0.1*[1 1 1];
            case 'legend'
                
                set(h.ItemText,'Color',[0.8 0.8 0.8]);
                h.TextColor = [0.8 0.8 0.8];
                h.Color = [0 0 0];
                
                if strcmp(fig.Tag,'SpectraPlot')
                    set(findall(fig,'Type','Patch'),'FaceAlpha',0.6);
                    set(findall(h.ItemTokens,'Type','Patch'),'FaceAlpha',0.6);
                end
                    
        end
    catch ME
        rethrow(ME)
    end
    
end


end