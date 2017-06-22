function goLight(fig)

fig.Color = [1 1 1];

childs = [findall(fig,'type','text');...
    findall(fig,'Style','text');...
    findall(fig,'type','axes');...
    findall(fig,'type','colorbar');...
    findall(fig,'type','uiflowcontainer');...
    findall(fig,'type','uicontainer'); ...
    findall(fig,'type','uigridcontainer'); ...
    findall(fig,'type','legend')];

for i = 1:length(childs)
    
    h = childs(i);
    
    if strcmp(h.Tag, 'ignore')
        continue;
    end
    
    try
        switch h.Style
            case 'text'
                h.ForegroundColor = [0.2 0.2 0.2];
                h.BackgroundColor = [1 1 1];
        end
    catch
    end
    
    try
        switch h.Type
            case 'text'
                
                h.Color = [0.2 0.2 0.2];
                
            case 'axes'
                
                h.Backdrop.FaceColor = [1 1 1];
                h.GridColor = [0.5 0.5 0.5];
                h.YColor = [0.2 0.2 0.2];
                h.XColor = [0.2 0.2 0.2];
                h.ZColor = [0.2 0.2 0.2];
                h.Title.Color = [0.1 0.1 0.1];
                h.Box = 'off';
                h.XRuler.TickDir = 'out';
                h.YRuler.TickDir = 'out';
                h.XGrid = 'on';
                h.YGrid = 'on';
                
            case 'colorbar'
                
                h.Ruler.Color = [0.2 0.2 0.2];
                h.BoxHandle.Visible = 'off';
                h.TickLength = 0.01;
                h.Ruler.TickDir = 'out';
                h.Ruler.Axle.Visible = 'off';
                
            case {'uiflowcontainer','uicontainer','uigridcontainer'}
                
                h.BackgroundColor = [1 1 1];
            case 'legend'
                set(h.ItemText,'Color',[0.2 0.2 0.2]);
                h.TextColor = [0.2 0.2 0.2];
                h.Color = [1 1 1];
                
                if strcmp(fig.Tag,'SpectraPlot')
                    set(findall(fig,'Type','Patch'),'FaceAlpha',0.3);
                    set(findall(h.ItemTokens,'Type','Patch'),'FaceAlpha',0.3);
                end
        end
    catch
    end
    
end



end