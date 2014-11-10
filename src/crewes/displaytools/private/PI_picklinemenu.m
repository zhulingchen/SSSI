function PI_picklinemenu();
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS

checkclick=get(gcf,'selectiontype');
hobj=gco;
if(strcmp(checkclick,'alt'))
    linenum=[];
    for ii=1:size(PICKS,1)
        CheckFigure=PICKS{ii};
        if(CheckFigure==gcf)
            CheckHandles=PICKS{ii,3};
            % fail safe
            if(isempty(CheckHandles))
                delete(hobj);
                return
            end
            for jj=1:size(CheckHandles)
                if(hobj==CheckHandles(jj))
                    linenum=jj;
                    break
                end
            end
        end
    end
    cmenu=uicontextmenu;
    set(gco,'uicontextmen',cmenu);
    % fail safe
    if(isempty(linenum))
        delete(hobj)
        return
    end
    PickMarker=findobj(gcf,'type','line','tag','PICKMARKER');
    m1=uimenu(cmenu,'label',['Move: ' num2str(linenum)],'callback','plotimage(''MovePickLine'')');
    m1=uimenu(cmenu,'label','Delete','callback','plotimage(''DeletePickLine'')');
elseif(strcmp(checkclick,'normal'))
    checkaction=get(gcf,'windowbuttondownfcn');
    if(strcmp(checkaction,'plotimage(''PickMoveClose'')'))
        % If the gcf propertie is above, pick lines are in move mode
    end
end
