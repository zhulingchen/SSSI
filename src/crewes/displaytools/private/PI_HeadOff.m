function StopOrGo=PI_HeadOff(action);
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
if(strcmp(get(mainax,'visible'),'off'))
    StopOrGo='STOP';
    return
end
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
h=get(gcf,'userdata');
hmsg=h(2);
CheckClick=get(gcf,'selectiontype');
StopOrGo='STOP';
if(strcmp(get(gca,'tag'),'POSITIONAXES'))
    box=selboxfini;
    if(length(box)==5)
        if(~ishandle(box(5)))
            return
        end
        delete(box(5));
    end
    return
end
if(strcmp(CheckClick,'normal'))
    if(strcmp(action,'PickMoveClose'))
        set(gcf,'name','Click & Hold MB1 on Markers to move.  MB3 menu, or function change to stop line moving');
        StopOrGo='STOP';
        return
    end
    StopOrGo='GO';
    % skipping
elseif(strcmp(get(gco,'type'),'image')&& strcmp(CheckClick,'alt'))
    h=get(gcf,'userdata');
    hzoompick=h(9);
    value=get(hzoompick,'value');
    delete(findobj(gcf,'type','line','tag','PICKMARKER'));
    delete(findobj(gcf,'type','text','tag','PICKTEXT'));
    switch value
    case 1
        selboxinit('plotimage(''zoom'')',1);
        set(gcf,'name','Seismic Image Plot, Simplezooming installed (Use MB1)');
    case 2
        drawlineinit('plotimage(''pick'')');
        set(gcf,'name','Seismic Image Plot, Picking resumed (Use MB1)');
    case 3
        drawlineinit('plotimage(''pick'')');
        set(gcf,'name','Seismic Image Plot, Picking new (Use MB1)');
        set(hzoompick,'userdata',[]);
    end
    set(gcf,'pointer','arrow');
    h=get(gcf,'userdata');
    hlimbox=h(14);
    val=get(hlimbox,'value');
    lmboxquestion={'Limit Box: ON' 'Limit Box: OFF'};
    sizequestion={'Publishable Figure: ON' 'Publishable Figure: OFF'};
    cmenu=uicontextmenu;
    set(gco,'UIContextMenu',cmenu);
    m1=uimenu(cmenu,'label','Zoom Options');
    m2=uimenu(m1,'label','Publish Zoom limits','callback','plotimage(''zoomoptions'')','userdata',[1 1]);
    jj=[1];
    PublishedFigure=[0];
    if(~isempty(ZOOM_VALUE))
        jj=[2];
        PublishedFigure=ZOOM_VALUE(5);
    end
    MenuLabel={'No Limits Published' 'Match Zoom To Published Limits'};
    en={'off' 'on'};
    m2=uimenu(m1,'label',MenuLabel{jj},'callback','plotimage(''zoomoptions'')','userdata',[1 2],...
        'enable',en{jj});
    CheckLockMat=[];
    MenuLabel={'Lock Zoom' 'Unlock Zoom'};
    UDat={[1 3] [1 4]};
    MenuCall={' ' 'plotimage(''zoomoptions'')'};
    jj=[1];
    for ii=1:size(ZOOM_LOCKS,1)
        CheckLock=ZOOM_LOCKS(ii,1);
        if(CheckLock==gcf)
            jj=[2];
            break
        end
    end
    m2=uimenu(m1,'label',MenuLabel{jj},'userdata',UDat{jj},'enable','on',...
        'callback',MenuCall{jj});
    if(jj==2)
    else
        PlotImageFigures=findobj(0,'type','figure','tag','PLOTIMAGEFIGURE');
        ifind=find(PlotImageFigures~=gcf);
        FigureLabel=PlotImageFigures(ifind);
        for ii=1:size(PICKS,1)
            if(PICKS{ii,1}==gcf)
            else
                haxs=findobj(PICKS{ii,1},'type','axes','tag','MAINAXES');
                ttl=get(haxs,'title');
                str=get(ttl,'string');
                m3=uimenu(m2,'label',['Figure: ' str],...
                    'callback','plotimage(''zoomoptions'')','userdata',[PICKS{ii} 3]);
            end
        end
    end
    ttl=get(mainax,'title');
    ttludat=get(ttl,'userdata');
    nm1=ttludat{3};
    m1=uimenu(cmenu,'label','Import Picks','enable','on');
    if(~isempty(PICKS))
        if(~isempty(PICKS{1,1}))
            set(m1,'enable','off');
            for ii=1:size(PICKS,1)
                if(PICKS{ii,1}==gcf)
                else
                    haxs=findobj(PICKS{ii,1},'type','axes','tag','MAINAXES');
                    ttl=get(haxs,'title');
                    ttludat=get(ttl,'userdata');
                    if(~iscell(ttludat))
                        
                    else
                        nm2=ttludat{3};
                        if(strcmp(nm1,nm2))
                            str=get(ttl,'string');
                            m2=uimenu(m1,'label',['Figure: ' str],...
                                'userdata',PICKS{ii,2},'callback','plotimage(''ImportPicks'')');
                            set(m1,'enable','on');
                        end
                    end
                end
            end
        else
            set(m1,'enable','off');
        end
    else
        set(m1,'enable','off');
    end
    m1=uimenu(cmenu,'label','Rename Axes','callback',@PI_axis_options,'separator','on');
    m1=uimenu(cmenu,'label','Resample Axes','callback',@PI_axis_options);
    
    m1=uimenu(cmenu,'label',lmboxquestion{val},'callback','plotimage(''LmLnActivation'')',...
        'separator','on');
    checkcolor=get(gcf,'color');
    if(checkcolor==[1 1 1])
        ii=2;
    else
        ii=1;
    end
    m1=uimenu(cmenu,'label','Send to Clipboard','callback',@PI_axis_options,'separator','on');
    if(strcmp(get(posax,'visible'),'on'))
        m1=uimenu(cmenu,'label','Hide Image Controls','callback',@PI_axis_options);
    else
        m1=uimenu(cmenu,'label','Show Image Controls','callback',@PI_axis_options);
    end
    m1=uimenu(cmenu,'label','Data Stats','callback',@PI_axis_options);
    m1=uimenu(cmenu,'label',sizequestion{ii},'callback','plotimage(''figuresizechange'')',...
        'separator','on');
    StopOrGo='STOP';
elseif(strcmp(get(gco,'type'),'image')&& strcmp(CheckClick,'extend'))
    % allowing for zoom motion
    h=get(gcf,'userdata');
    hi=h(5);
    dat=get(hi,'cdata');
    if(size(dat,1)*size(dat,2)>=600000)
        StopOrGo='STOP';
        return
    end
    plotimage('zoomscroll');
    StopOrGo='STOP';
elseif(strcmp(CheckClick,'alt'))
    % building uicontext menu for Spawned Plot Image
    cmenu=uicontextmenu;
    set(gcf,'UIContextMenu',cmenu);
    sep='on';
    if(length(findobj(gcf,'type','uimenu','tag','PLOTIMAGEMASTERFIGURE'))==1)
        sep='on';
        m1=uimenu(cmenu,'label','Change Globals','callback',@PI_Global);
        m1=uimenu(cmenu,'label','Spawn New Plot Image','callback','plotimage(''SpawnPlotImage'')','separator',sep);
    end
    checkcolor=get(gcf,'color');
    if(checkcolor==[1 1 1])
        ii=2;
    else
        ii=1;
    end
    m1=uimenu(cmenu,'label','Send to Clipboard','callback',@PI_axis_options,'separator','on');
    if(strcmp(get(posax,'visible'),'on'))
        m1=uimenu(cmenu,'label','Hide Image Controls','callback',@PI_axis_options);
    else
        m1=uimenu(cmenu,'label','Show Image Controls','callback',@PI_axis_options);
    end
    m1=uimenu(cmenu,'label','Data Stats','callback',@PI_axis_options);
    sizequestion={'Publishable Figure: ON' 'Publishable Figure: OFF'};
    m1=uimenu(cmenu,'label',sizequestion{ii},'callback','plotimage(''figuresizechange'')',...
        'separator',sep);
    
end
