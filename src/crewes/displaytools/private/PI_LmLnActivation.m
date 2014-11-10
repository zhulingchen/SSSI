function PI_LmLnActivation
% -------------------------------
% ----- Limit Box Callbacks -----
% -------------------------------
%
% This ssction deals with the call backs for the Limit box that acts
% simialr to the slave and mast routin that is present plotimage
% A hidden button is being used to store data for the lines
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
h=get(gcf,'userdata');
hmsg=h(2);
hi=h(5);
hscale=h(6);
hclip=h(7);
hmaster=h(10);
hlimbox=h(14);
ckval=get(hmaster,'value');
cklimbox=get(hlimbox,'value');
if(cklimbox==1)
    set(hlimbox,'value',[2]);
else
    set(hlimbox,'value',[1])
end
if(ckval==3)
    % figure has been slaved to another, asking user to unslave or not
    checkwuser=questdlg('Do you want to unslave present figure?',...
        'Slave Alert','Yes','No','Cancel','Cancel');
    switch checkwuser
    case 'Yes'
        set(hmaster,'value',[1]);
    case 'No'
        set(hlimbox,'value',[1]);
        return   
    case 'Cancel'
        set(hlimbox,'value',[1]);
        return
    end
end
chval=get(hlimbox,'value');
% 1 - Limit box will be shut off
% 2 - Limit Box will be activated
limdat=get(hlimbox,'userdata');
checklns=findobj(gcf,'type','line','tag','LIMITLINE');
if(~isempty(limdat)&~isempty(checklns))
    % setting limit box off
    if(chval==1)
        set(limdat{1},'visible','off');
        set(limdat{2},'visible','off');
        set(limdat{4},'visible','off');
        set(findobj(gcf,'type','uimenu','tag','LIMITBOXMENU'),'label','Limit Box: ON');
        stringinfo='Visibility of Limits lines has been shut off';
    elseif(chval==2)
        set(findobj(gcf,'type','uimenu','tag','LIMITBOXMENU'),'label','Limit Box: OFF');
        set(limdat{1},'visible','on');
        set(limdat{2},'visible','on');
        set(limdat{4},'visible','on');
        stringinfo='Visibility of Limits lines has been shut on';
    end
else
    % need to build lines, points and measurementfig
    axes1=mainax;
    xdat=get(axes1,'xlim');
    xmin=xdat(1);
    xmax=xdat(2);
    ydat=get(axes1,'ylim');
    ymin=ydat(1);
    ymax=ydat(2);
    col=[1 0 1];
    limlnstyle='--';
    limlnbd='plotimage(''limlnmove'')';
    limlntag='LIMITLINE';
    limptmk='.';
    limptmksz=[18];
    limptbd='plotimage(''limptmove'')';
    limpttag='LIMITPOINT';
    lwid=[2];
    ddy=get(hi,'ydata');
    ddx=get(hi,'xdata');
    xmin1=xmin; xmax1=xmax; ymin1=ymin; ymax1=ymax;
    dy=(ydat(2)-ydat(1))/120;
    ymin=max([ymin+dy ddy(1)]);
    ymax=min([ymax-dy ddy(end)]);
    dx=(xdat(2)-xdat(1))/120;
    xmin=max([xmin+dx ddx(1)]);
    xmax=min([xmax-dx ddx(end)]);
    dy=0;
    dx=0;
    % plotting lines
    % bottom
    limln1=line([xmin xmax],[ymax-dy ymax-dy],'linestyle',limlnstyle,'buttondownfcn',...
        limlnbd,'color',col,'tag',limlntag,'linewidth',lwid);
    % top
    limln2=line([xmin xmax],[ymin+dy ymin+dy],'linestyle',limlnstyle,'buttondownfcn',...
        limlnbd,'color',col,'tag',limlntag,'linewidth',lwid);
    % Left Side
    limln3=line([xmin+dx xmin+dx],[ymin ymax],'linestyle',limlnstyle,'buttondownfcn',...
        limlnbd,'color',col,'tag',limlntag,'linewidth',lwid);
    % Right Side
    limln4=line([xmax-dx xmax-dx],[ymin ymax],'linestyle',limlnstyle,'buttondownfcn',...
        limlnbd,'color',col,'tag',limlntag,'linewidth',lwid);
    
    % Plotting Points
    % Upper Left
    limpt1=line([xmin+dx],[ymax-dy],'marker',limptmk,'buttondownfcn',limptbd,'color',col,'tag',limpttag,...
        'markersize',limptmksz);
    % Upper Right
    limpt2=line([xmax-dx],[ymax-dy],'marker',limptmk,'buttondownfcn',limptbd,'color',col,'tag',limpttag,...
        'markersize',limptmksz);
    % Lower left
    limpt3=line([xmin+dx],[ymin+dy],'marker',limptmk,'buttondownfcn',limptbd,'color',col,'tag',limpttag,...
        'markersize',limptmksz);
    % Lower Right
    limpt4=line([xmax-dx],[ymin+dy],'marker',limptmk,'buttondownfcn',limptbd,'color',col,'tag',limpttag,...
        'markersize',limptmksz);
    
    % Center Point
    limcent=line((xmax1-xmin1)/2+xmin1,(ymax1-ymin1)/2+ymin1,'marker','x','buttondownfcn','plotimage(''limcentmove'')','color',col,...
        'tag',limpttag,'userdata',[3],'markersize',[8]);
    
    % [upper left upper right lower left lower right] [top bottom left side right side]
    limlnpositions=[limln1 limln2 limln3 limln4];
    limptpositions=[limpt1 limpt2 limpt3 limpt4];
    hmasterfig=gcf;
    mfig=findobj(0,'type','figure','tag','MVLINESMEASUREMENTS');
    buildfigure=1;
    for ii=1:length(mfig)
        checkfig=get(mfig(ii),'userdata');
        if(checkfig==hmasterfig)
            buildfigure=[];
            break
        end
    end
    if(~isempty(buildfigure))
        measurementfig=figcent(.1,.21);
        set(measurementfig,'menubar','none','tag','MVLINESMEASUREMENTS',...
            'numbertitle','off','name',['Limit for fig: ' num2str(hmasterfig) ],'buttondownfcn',' ',...
            'visible','off','closerequestfcn','plotimage(''limlnfigurevis2'')',...
            'resize','off','userdata',hmasterfig);
        nnam=uicontrol(measurementfig,'style','Text','string','Limit Lines Data','units','normalized',...
            'position',[0 .858 1 .13],'tag','MEASURENAME','backgroundcolor',[0 0 0],...
            'foregroundcolor',[1 1 1]);
        xx=uicontrol(measurementfig,'style','Text','string','Bottom','units','normalized',...
            'position',[-.007 0.747 0.340 .108],'tag','MEASURETOP');
        xx=uicontrol(measurementfig,'style','Text','string','Top','units','normalized',...
            'position',[-.007 0.626 0.340 .108],'tag','MEASUREBOTTOM');
        xx=uicontrol(measurementfig,'style','Text','string','Left','units','normalized',...
            'position',[-.007 0.506 0.340 .108],'tag','MEASURELEFT');
        xx=uicontrol(measurementfig,'style','Text','string','Right','units','normalized',...
            'position',[-.007 0.386 0.340 .108],'tag','MEASURERIGHT');
        mbot=uicontrol(measurementfig,'style','Edit','string',num2str(ymin+dy),'units','normalized',...
            'position',[0.347 0.753 0.64 .096],'tag','MEASURETOPNUM',...
            'backgroundcolor',[1 1 1]);
        mtop=uicontrol(measurementfig,'style','Edit','string',num2str(ymax-dy),'units','normalized',...
            'position',[0.347 0.633 0.64 .096],'tag','MEASUREBOTTOMNUM',...
            'backgroundcolor',[1 1 1]);
        mlft=uicontrol(measurementfig,'style','Edit','string',num2str(xmin+dx),'units','normalized',...
            'position',[0.347 0.512 0.64 .096],'tag','MEASURELEFTNUM',...
            'backgroundcolor',[1 1 1]);
        mrit=uicontrol(measurementfig,'style','Edit','string',num2str(xmax-dx),'units','normalized',...
            'position',[0.347 0.392 0.64 .096],'tag','MEASURERIGHTNUM',...
            'backgroundcolor',[1 1 1]);
        xx=uicontrol(measurementfig,'string','Reset','units','normalized',...
            'position',[0 .229 .340 .139],'tag','MEASURERESET','callback','plotimage(''lmptreset'')');
        xx=uicontrol(measurementfig,'string','Enter','units','normalized',...
            'position',[.353 .229 .633 .139],'tag','MEASUREENTER','callback','plotimage(''lmptenter'')');
        mhdf=uicontrol(measurementfig,'style','togglebutton','string','Hide Data Figure',...
            'units','normalized','position',[0.019 0 0.962 0.13],'tag','MEASUREHIDE',...
            'callback','plotimage(''limlnfigurevis2'')','max',[1],'min',[2],'value',[1]); 
        % saving handles
        mh=[mbot mtop mlft mrit mhdf hmasterfig];
        set(measurementfig,'userdata',mh);
        limdata={limptpositions limlnpositions measurementfig limcent};
        set(hlimbox,'userdata',limdata);
        stringinfo='Limit lines have been actiavted';
    else
    end
    set(findobj(gcf,'type','uimenu','tag','LIMITBOXMENU'),'label','Limit Box: OFF');
end
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
