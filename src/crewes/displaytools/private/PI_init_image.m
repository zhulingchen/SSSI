function PI_init_image()
global SCALE_OPT AMPFLAG NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP CLOSEREQUEST
global CHECK_GLOBAL
if(isempty(CHECK_GLOBAL))
    % user has not predefined 
    try
        load('plotimageproperties.mat');
        sv=2;
    catch
        sv=1;
    end
else
    % user has predefined globals
    sv=1;
end
if(isempty(AMPFLAG)); Ampflag='I'; else Ampflag=AMPFLAG; end
if(isempty(NOBRIGHTEN)); nobrighten=0; else nobrighten=NOBRIGHTEN; end
if (isempty(SCALE_OPT)); scaleopt=2; else scaleopt=SCALE_OPT; end %default to max scaling
if (isempty(NUMBER_OF_COLORS)); number_of_colors=64;
else number_of_colors=NUMBER_OF_COLORS; end %default to 64 gray levels
if (isempty(GRAY_PCT)); gray_pct=50; else gray_pct=GRAY_PCT; end %default to 50% gray transition
if (isempty(CLIP)); clip=4; else clip=CLIP; end %default to clip of 4
if (isempty(COLOR_MAP)|strcmp(COLOR_MAP,'seisclrs'))
    clrmap=seisclrs(number_of_colors,gray_pct);
else
    clrmap=COLOR_MAP;
end
if (isempty(PICKCOLOR))
    PICKCOLOR='r';
end
if (isempty(XAXISTOP))
    XAXISTOP=0;
end
if(isempty(CLOSEREQUEST))
    CLOSEREQUEST='Slow Close';
end
    
if(sv==1)
    save('plotimageproperties.mat','XAXISTOP','PICKCOLOR','NOBRIGHTEN','COLOR_MAP',...
    'CLIP','GRAY_PCT','NUMBER_OF_COLORS','SCALE_OPT','AMPFLAG');
end
gray_pct=round(gray_pct); %force an integerscaleopt=[];
mxs2=[]; % set(hscale,'userdata',[scaleopt mxs2 mns smean2 stddev2]);
mns=[]; % set(hmaster,'userdata',[mxs smean stddev]);
smean2=[];
stddev2=[];
mxs=[];
smean=[];
stddev=[];
hi=[0];     % no image yet
switch Ampflag
    case 'I'
        ampflag=1;
    case 'M'
        ampflag=2;
    case 'S'
        ampflag=3;
end
clips=[30 25 20 15 10 9 8 7 6 5 4 3 2 1 .5 .25 .1 .05 .01 .005 .001];
iclip=near(clips,clip);
clip=clips(iclip);
enb='off';   % buttons not allowed on yet
% Building menus
checkscreen=get(0,'screensize');
if(checkscreen(3)>=1000)
    wid=.4;
else
    wid=.8;
end
cfig1=figcent(wid,.5);
thisfig=gcf;
set(cfig1,'name',['Seismic Image Plot(' int2str(thisfig) '), Simplezooming installed (Use MB1)'],...
   'closerequestfcn',@PI_Close,'tag','PLOTIMAGEFIGURE',...
   'menubar','none','numbertitle','off');
hposax=subplot('position',[.029 .721 .158 .18]);
set(hposax,'hittest','off','visible','off','tag','POSITIONAXES');
hgca=subplot('position',[.291 .117 .686 .783]);
set(hgca,'visible','off','tag','MAINAXES');
selboxinit('plotimage(''zoom'')',1);
% File Menus
hfile1=uimenu(gcf,'Label','File','enable','on');
hf=uimenu(hfile1,'Label','Open .mat','visible','on','callback','plotimage(''OpenFile'')',...
    'userdata',{[1] gcf []});
hf=uimenu(hfile1,'Label','Open Segy','visible','on','callback','plotimage(''OpenFile'')',...
    'userdata',{[2] gcf []});
hfsave=uimenu(hfile1,'label','Save','callback','plotimage(''SaveFile'')','separator','on',...
    'userdata',{[1] gcf []});
hf=uimenu(hfile1,'label','Save As','callback','plotimage(''SaveFile'')',...
    'userdata',{[2] gcf []});
% hf=uimenu(hfile1,'label','Open Properties','callback','ChangePropertiesMenu',...
%     'userdata',{[4] gcf []},'separator','on');
% hf=uimenu(hfile1,'label','Save Properties','callback','plotimage(''SaveFile'')',...
%     'userdata',{[3] gcf []});
if(length(findobj(0,'type','figure','tag','PLOTIMAGEFIGURE'))<=1)
    hf=uimenu(hfile1,'label','Spawn New Plot Image','callback','plotimage(''SpawnPlotImage'')','separator','on',...
        'tag','PLOTIMAGEMASTERFIGURE');
else
    masterfigchild=findobj(0,'type','uimenu','tag','PLOTIMAGEMASTERFIGURE');
    if(isempty(masterfigchild))
        hf=uimenu(hfile1,'label','Spawn New Plot Image','callback','plotimage(''SpawnPlotImage'')','separator','on',...
        'tag','PLOTIMAGEMASTERFIGURE');
    else
        masterfigchild=get(masterfigchild,'parent');
        masterfig=get(masterfigchild,'parent');
        pos=get(masterfig,'position');
        set(gcf,'position',[pos(1)+20 pos(2)-20 pos(3) pos(4)]);
    end    
end
filenames=[];
lastpath=pwd;
originalpath=pwd;
try
    load('plotimagedata.mat');
catch
end
previousload=[];
if(isempty(filenames))
    % if there is no previous loaded session, this is creating blank
    % menu items
    xx1=uimenu(hfile1,'separator','on','callback','plotimage(''OpenFile'')',...
        'visible','off','userdata',{[3] gcf []},'tag','QUICK_OPENFILE');
    previousload=[previousload xx1];
    for ii=1:3
        xx=uimenu(hfile1,'callback','plotimage(''OpenFile'')',...
            'visible','off','userdata',{[3] gcf []});
        previousload=[previousload xx];
    end
    set(xx1,'userdata',{[3] gcf previousload originalpath lastpath});
elseif(~isempty(filenames))
    % this is creating menu items for the files that were previously load
    % in other sessions
    for ii=1:4
        vis='on';
        if(isempty(strunpad(filenames{ii,:})))
            vis='off';
        end
        if(size(filenames,1)>=ii)
            lbl=strunpad(filenames{ii,:});
            xx1=uimenu(hfile1,'label',strunpad(lbl),...
                'callback','plotimage(''OpenFile'')','visible',vis,...
                'userdata',{[3] gcf []});
        else
            xx1=uimenu(hfile1,'callback','plotimage(''OpenFile'')',...
                'visible',vis,'userdata',{[3] gcf []},'userdata');
        end
        previousload=[previousload xx1];
        set(previousload(1),'userdata',{[3] gcf previousload originalpath lastpath},...
            'separator','on','tag','QUICK_OPENFILE');
    end
end
hf=uimenu(hfile1,'label','Close','callback',@PI_Close,'separator','on',...
    'userdata',[gcf]);  % closing will be different for master figs and spawened figs
% Option Menus
hfile2=uimenu(gcf,'label','Options','tag','PLOTIMAGEOPTIONS','enable',enb);
hf1=uimenu(hfile2,'label','Limit Box: ON','enable','on','tag','LIMITBOXMENU',...'
    'callback','plotimage(''LmLnActivation'')');
hf2=uimenu(hfile2,'label','Visibility of Data Figure','callback','plotimage(''limlnfigurevis'')',...
    'enable','off','tag','LIMITBOXDATAFIGUREMENU','visible','off');
if(length(findobj(0,'type','figure','tag','PLOTIMAGEFIGURE'))<=1)
%     hf3=uimenu(hfile2,'label','Change Properties','callback','plotimage(''ChangeProperties'')',...
%         'separator','on');
    
    mrk=uimenu(gcf,'label',['              Plot Image Master Figure ' date]);
end
hf=uimenu(hfile2,'label','Open Pick File','callback','plotimage(''PicksOpen'')','separator','on');
hf=uimenu(hfile2,'label','Save Picks','callback','plotimage(''PicksSave'')');
hamp_picks=uimenu(hfile2,'label','Copy AMP_PICKS from another window','callback','plotimage(''copyamppicks'')');
hamp_picks_lbl=uimenu(hfile2,'label','Show AMP_PICKS event names','callback','plotimage(''showamppicks'')','checked','on');
    %%%% Margrave May 30 2006 Hack begins
    posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
    sep='on';
    if(length(findobj(gcf,'type','uimenu','tag','PLOTIMAGEMASTERFIGURE'))==1)
        sep='on';
        m1=uimenu(hfile2,'label','Change Globals','separator','on','callback',@PI_Global);
        m1=uimenu(hfile2,'label','Spawn New Plot Image','callback','plotimage(''SpawnPlotImage'')','separator',sep);
    end
    checkcolor=get(gcf,'color');
    if(checkcolor==[1 1 1])
        ii=2;
    else
        ii=1;
    end
    m1=uimenu(hfile2,'label','Send to Clipboard','callback',@PI_axis_options,'separator','on');
    if(strcmp(get(posax,'visible'),'on'))
        m1=uimenu(hfile2,'label','Hide Image Controls','callback',@PI_axis_options);
    else
        m1=uimenu(hfile2,'label','Hide Image Controls','callback',@PI_axis_options);
    end
    m1=uimenu(hfile2,'label','Data Stats','callback',@PI_axis_options);
    sizequestion={'Publishable Figure: ON' 'Publishable Figure: OFF'};
    m1=uimenu(hfile2,'label',sizequestion{ii},'callback','plotimage(''figuresizechange'')',...
        'separator',sep);
    %%%%% end of Margrave Hack

if(isempty(PICKS))
    PICKS={[gcf] [] []};
else
    PICKS{size(PICKS,1)+1,1}=gcf;
end

%put the x axis on top
if(XAXISTOP)
    set(gca,'xaxislocation','top')
end

%make a few buttons
enb='off';
% Text backing
hbak=uicontrol('style','text','string','Image Controls',...
    'units','normalized','position',[.018 .117 .18 .588],...
    'backgroundcolor',[.502 .502 .502],'foregroundcolor',[1 1 1],...
    'tag','BACKING');

hzoompick=uicontrol('style','popupmenu',...
    'string',['Zoom|Pick time dips(Old PICK buffer)|Pick time dips(New PICK buffer)'...
        '|Pick amplitudes(Old AMP_PICK buffer)|Pick amplitudes(New AMP_PICK buffer)'],...
    'units','normalized','tooltipstring','Define mouse action as zoom or pick',...
    'position',[.018 .605 .18 .05],'callback','plotimage(''zoompick'')',...
    'enable',enb);
hflip=uicontrol('style','popupmenu',...
    'string','Normal Polarity|Reverse Polarity',...
    'units','normalized','tooltipstring','Set display polarity',...
    'position',[.018 .543 .18 .05],'callback','plotimage(''flip'')',...
    'userdata',1,'enable',enb);
fsize=get(0,'factoryuicontrolfontsize');
hlabel=uicontrol('style','text','fontsize',fsize,'string','Bright 0','units','normalized',...
    'position', [.018 .488 .18 .033],'tooltipstring','Current brightness level','userdata',0,'enable',enb);
hslider=uicontrol('style','slider','string','Bright','units','normalized','position',...
    [.018 .463 .18 .025],'callback','plotimage(''brighten'')',...
    'tooltipstring','Set image brightness','max',10,'min',-10,...
    'tag','phase','value',0,'userdata',hlabel,'sliderstep',[1/20 1/20],'enable',enb);
% user data for below is being used
hmsg=uicontrol('style','text','string','Polarity Normal',...
    'units','normalized',...
    'position',[0 0 .15 .05],'visible','off','enable',enb);
hmaster=uicontrol('style','popupmenu','string','Independent|Master|Slave',...
    'units','normalized','position',[.018 .398 .18 .05],'tooltipstring','Define amplitude control',...
    'callback','plotimage(''rescale'')','value',ampflag,'enable',enb);
hscale=uicontrol('style','popupmenu','string',str2mat('Mean scaling',...
    'Max scaling'),'units','normalized','position',[.018 .336 .18 .05],'tooltipstring','Define data scaling mechanism',...
    'callback','plotimage(''rescale'')','value',scaleopt,'enable',enb);
vis='on'; if(scaleopt==2) vis='off'; end
nclips=length(clips);
clipmsg=num2strmat(clips);
clipmsg=[ones(nclips,1)*'Cliplevel: ' num2strmat(clips)];
hclip=uicontrol('style','popupmenu','string',clipmsg,...
    'units','normalized','position',[.018 .273 .18 .05],'tooltipstring','Set clip level in std deviations',...
    'callback','plotimage(''rescale'')','value',iclip,...
    'visible',vis,'enable',enb,'userdata',clips);
boxstr=strmat('Limit Box: OFF','Limit Box: ON ');
hlimbox=uicontrol('style','popupmenu','string',boxstr,...
    'units','normalized','position',[0 0 .15 .05],...
    'tooltipstring','Constrains master values to box',...
    'callback','LmLnActivation','value',[1],...
    'visible','off','enable',enb);
% settign all buttons into for quick access
set(hbak,'userdata',[hbak hzoompick hflip hlabel hslider hmaster hscale hclip]);
% hlimbox visibility is now set off and actions will be accessed
% through a uicontext menu.  This context menu is built in PI_HeadOff. THis
% is a bit of a problem since the context menu is an attribute of the image
% and the image is supposed to respond to clicks for zooming and picking.
% Since MB3 invokes the context menu but is also involved in picking, this
% is a bit of a conflict. At present, the context menu only displays in
% "zoom" mode.

% zoom controls
hvrtscrol=uicontrol('style','slider','units','normalized','position',[.957 .139 .021 .762],...
    'callback',@PI_zoom_slider,'userdata',{[1] [] []},'visible','off','sliderstep',[.05 .2]);
hhorscrol=uicontrol('style','slider','units','normalized','position',[.291 .117 .664 .021],...
    'callback',@PI_zoom_slider,'userdata',{[2] [] []},'visible','off','sliderstep',[.05 .2]);

% message box, 
hmsgmain=uicontrol('style','text','units','normalized','string','Welcome to Plotimage',...
    'position',[0 0 1 .035],'backgroundcolor',[1 1 1],'visible','on');

% position axes
hposax=subplot('position',[.029 .721 .158 .18]);
set(hposax,'hittest','off','visible','off','tag','POSITIONAXES');
set(gcf,'currentaxes',hgca);
 
set(gcf,'userdata',[hflip hmsgmain -1 hmsg hi hscale hclip -1 ,...
        hzoompick hmaster hlabel hslider hclip hlimbox hfile2 hvrtscrol hhorscrol hamp_picks_lbl]);
    
    set(hscale,'userdata',[scaleopt mxs2 mns smean2 stddev2]);
    set(hmaster,'userdata',[mxs smean stddev]);
    set(hclip,'userdata',iclip);
    set(hmsg,'userdata',clips)
    
global IMCONTROLS
% this was add on Ocotber 27th for Gary Margrave to allow users to
% automatically shup off the gui controls.
if(isempty(IMCONTROLS))
    IMCONTROLS='on';
end
if(strcmp(lower(IMCONTROLS),'on'))
    % do nothing
elseif(strcmp(lower(IMCONTROLS),'off'))
    PI_axis_options;
end
