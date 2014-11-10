function PI_limmoveend();
% this routine occures at the end of moving the limit lines.  It ensure
% that the limit box parts are reset to the original order.  This section
% also forces the limit lines to be 5% away from the edge of the image
% data.
global PICKS
if(~isempty(PICKS))
    for ii=1:size(PICKS,1)
        CheckFigure=PICKS{ii,1};
        if(CheckFigure==gcf)
            set(PICKS{ii,3},'erasemode','normal');
        end
    end
end
h=get(gcf,'userdata');
hi=h(5); % Image
hlimbox=h(14);
limdata=get(hlimbox,'userdata');
% setting points and lines erase mode to normal
limlns=limdata{2};
limpts=limdata{1};
limcent=limdata{4};
set(limlns,'erasemode','normal');
set(limpts,'erasemode','normal','visible','on');
set(limcent,'erasemode','normal');
set(gcf,'windowbuttonupfcn',[]);
set(gcf,'windowbuttonmotionfcn',[]);
selboxinit('plotimage(''zoom'')',1);
h=get(gcf,'userdata');
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
limptdat=limdat{1};
limlndat=limdat{2};
findfig=limdat{3};
% need to stop axis from resizing
parentaxis=get(limptdat(1),'parent');
xdat=get(parentaxis,'xlim');
ydat=get(parentaxis,'ylim');
% [top bottom left right]
ximdat=get(hi,'xdata');
yimdat=get(hi,'ydata');
tpnum=get(limlndat(1),'ydata');
btnum=get(limlndat(2),'ydata');
dy=sort([ydat(1) ydat(2)]);
dy=(dy(2)-dy(1))/120;
hornum=sort([yimdat(1)+dy tpnum(1) btnum(2) yimdat(end)-dy]);
hornum=[hornum(2) hornum(3)];
lfnum=get(limlndat(3),'xdata');
rtnum=get(limlndat(4),'xdata');
dx=sort([xdat(1) xdat(2)]);
dx=(dx(2)-dx(1))/120;
vernum=sort([ximdat(1)+dx lfnum(1) rtnum(1) ximdat(end)-dx]);
vernum=[vernum(2) vernum(3)];
mdat=[hornum vernum];   
mlninfo=get(findfig,'userdata');
ddat={'ydata' 'ydata' 'xdata' 'xdata'};
for ii=1:4
    set(mlninfo(ii),'string',num2str(mdat(ii)));
end
checkvis=mlninfo(5);
% setting lines back
set(limlndat(1),'ydata',[mdat(2) mdat(2)],'xdata',[mdat(3) mdat(4)]);
set(limlndat(2),'ydata',[mdat(1) mdat(1)],'xdata',[mdat(3) mdat(4)]);
set(limlndat(3),'ydata',[mdat(1) mdat(2)],'xdata',[mdat(3) mdat(3)]);
set(limlndat(4),'ydata',[mdat(1) mdat(2)],'xdata',[mdat(4) mdat(4)]);
% setting points back
set(limptdat(1),'ydata',mdat(2),'xdata',mdat(3));
set(limptdat(2),'ydata',mdat(2),'xdata',mdat(4));
set(limptdat(3),'ydata',mdat(1),'xdata',mdat(3));
set(limptdat(4),'ydata',mdat(1),'xdata',mdat(4));
set(parentaxis,'ylim',ydat,'xlim',xdat);
% setting center line back (if it needs to go back)
set(limcent,'ydata',(mdat(2)-mdat(1))/2+mdat(1),'xdata',(mdat(4)-mdat(3))/2+mdat(3));
if(get(checkvis,'value')==2)
    % measurement data figure can be made visible
    set(findfig,'visible','on');
else
    % measurement data figure will not be made visible
end
hgcf=gcf;
plotimage('limboxrescale');
set(0,'currentfigure',hgcf);
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
    set(gcf,'name','Seismic Image Plot, Picking resummed (Use MB1)');
case 3
    drawlineinit('plotimage(''pick'')');
    set(gcf,'name','Seismic Image Plot, Picking new (Use MB1)');
    set(hzoompick,'userdata',[]);
end
