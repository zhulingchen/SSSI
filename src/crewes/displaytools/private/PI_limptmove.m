function PI_limptmove(action);
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS

delete(findobj(gcf,'type','line','tag','PICKMARKER'));
h=get(gcf,'userdata');
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
st=get(gcf,'selectiontype');
if(strcmp(st,'alt'))
    % right click, creating menu
    axes1=findobj(gcf,'type','axes','tag','AXES1');
    if(ishandle(get(axes1,'userdata')));
        delete(get(axes1,'userdata'));
    end
    cmenu=uicontextmenu;
    set(gco,'uicontextmen',cmenu);
    m1=uimenu(cmenu,'label','Line Options','callback','plotimage(''mvmenuclose'')');
    m2=uimenu(cmenu,'label','Line Style','separator','on');
    m2x=uimenu(m2,'label','- (Solid line)','callback','plotimage(''limlnoptions'')','userdata',[1 1]);
    m2x=uimenu(m2,'label','-- (Segmented line)','callback','plotimage(''limlnoptions'')','userdata',[1 2]);
    m2x=uimenu(m2,'label','_. (Line dot)','callback','plotimage(''limlnoptions'')','userdata',[1 3]);
    m2x=uimenu(m2,'label',': (Dotted line)','callback','plotimage(''limlnoptions'')','userdata',[1 4]);
    m3=uimenu(cmenu,'label','Marker');
    m3x=uimenu(m3,'label','*','callback','plotimage(''limlnoptions'')','userdata',[2 1]);
    m3x=uimenu(m3,'label','x','callback','plotimage(''limlnoptions'')','userdata',[2 2]);
    m3x=uimenu(m3,'label','o','callback','plotimage(''limlnoptions'')','userdata',[2 3]);
    m3x=uimenu(m3,'label','+','callback','plotimage(''limlnoptions'')','userdata',[2 4]);
    m3x=uimenu(m3,'label','.','callback','plotimage(''limlnoptions'')','userdata',[2 5]);
    m4=uimenu(cmenu,'label','Color');
    m4x=uimenu(m4,'label','Red','callback','plotimage(''limlnoptions'')','userdata',[3 1]);
    m4x=uimenu(m4,'label','Green','callback','plotimage(''limlnoptions'')','userdata',[3 2]);
    m4x=uimenu(m4,'label','Blue','callback','plotimage(''limlnoptions'')','userdata',[3 3]);
    m4x=uimenu(m4,'label','Black','callback','plotimage(''limlnoptions'')','userdata',[3 4]);
    m4x=uimenu(m4,'label','Default','callback','plotimage(''limlnoptions'')','userdata',[3 5]);
    m5=uimenu(cmenu,'label','Spawn Options','separator','on');
    ttl=get(gca,'title');
    ttludat=get(ttl,'userdata');
    nm1=ttludat{3};
    if(strcmp(nm1,'Amp Spectrum'))
        m5x=uimenu(m5,'label','None Yet For Spectrums');
    else
        m5x=uimenu(m5,'label','Spawn New Figure','callback',@PI_SpawnPlotImage,...
            'userdata',[]);
        m5x=uimenu(m5,'label','F-K Amp Spectrum','callback',@PI_limbox_transforms,...
            'userdata',[],'separator','on');
        m5x=uimenu(m5,'label','F-X Amp Spectrum','callback',@PI_limbox_transforms,...
            'userdata',[]);
    end
    m6=uimenu(cmenu,'label','Reset Data','callback','plotimage(''lmptresetmenu'')');
    slavefig=limdat{3};
    checkvis=get(slavefig,'visible');
    if(strcmp(checkvis,'on'))
        nn=2;
    else
        nn=1;
    end
    menstr={'Data Figure: ON' 'Data Figure: OFF'};
    m7=uimenu(cmenu,'label',menstr{nn},'callback','plotimage(''limnoptions'')',...
        'separator','on','callback','plotimage(''limlnfigurevis'')');
    set(axes1,'userdata',cmenu);
else
    if(~isempty(PICKS))
        for ii=1:size(PICKS,1)
            CheckFigure=PICKS{ii,1};
            if(CheckFigure==gcf)
                checkpics=findobj(gcf,'type','line','tag','PICKS');
                if(isempty(checkpics))
                    PICKS{ii,2}=[];
                    PICKS{ii,3}=[];
                    break
                end
                set(PICKS{ii,3},'erasemode','xor');
            end
        end
    end
    axes1=gca;
    % setting points and lines erase mode to xor
    lmlm=findobj(gcf,'type','uicontrol','tag','LIMITLINEMASTER');
    h=get(gcf,'userdata');
    hlimbox=h(14);
    limdat=get(hlimbox,'userdata');
    limpts=limdat{1};
    limlns=limdat{2};
    limcent=limdat{4};
    set(limlns,'erasemode','xor');
    set(limpts,'erasemode','xor');
    set(limcent,'erasemode','xor');
    set(gcf,'windowbuttonupfcn','plotimage(''limmoveend'')');
    if(strcmp(action,'limptmove'))
        set(gcf,'windowbuttonmotionfcn','plotimage(''limptmove2'')');
        % need to store line positions [top bottom left side right side]
        pts=cell(5,1);
        for ii=1:4
            xdat=get(limlns(ii),'xdata');
            ydat=get(limlns(ii),'ydata');
            tdat=[xdat ydat];
            pts{ii}=tdat;
        end
        pts{5}=get(gca,'currentpoint');
        set(axes1,'userdata',pts);
    elseif(strcmp(action,'limlnmove'))
        set(gcf,'windowbuttonmotionfcn','plotimage(''limlnmove2'')');
        pt=get(axes1,'currentpoint');
        set(axes1,'userdata',pt);
    elseif(strcmp(action,'limcentmove'))
        set(gcf,'windowbuttonmotionfcn',@PI_limcentmove2);
        pt=get(axes1,'currentpoint');
        set(axes1,'userdata',pt);
    end
    set(gcf,'windowbuttonupfcn','plotimage(''limmoveend'')');
end
