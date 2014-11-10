function PI_positionaxes_lineposition(hObject, eventdata, handles)
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
h=get(get(posax,'title'),'userdata');
% [bottom top left right]
xdat=get(mainax,'xlim');
ydat=get(mainax,'ylim');
if(strcmp(get(posax,'visible'),'on'))
    vis='on';
else
    vis='off';
end
set(h(1),'ydata',[ydat(1) ydat(1)],'xdata',[xdat(1) xdat(2)],'visible',vis);
set(h(2),'ydata',[ydat(2) ydat(2)],'xdata',[xdat(1) xdat(2)],'visible',vis);
set(h(3),'ydata',[ydat(1) ydat(2)],'xdata',[xdat(1) xdat(1)],'visible',vis);
set(h(4),'ydata',[ydat(1) ydat(2)],'xdata',[xdat(2) xdat(2)],'visible',vis);
set(h(6),'ydata',[ydat(1) ydat(2) ydat(2) ydat(1)],'xdata',[xdat(2) xdat(2) xdat(1) xdat(1)],'visible',vis);
