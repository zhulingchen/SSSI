function PI_limcentmove2(hObject, eventdata, handles)
axes1=gca;
h=get(gcf,'userdata');
hi=h(5);
hlimbox=h(14);

xlimits=get(hi,'xdata');
xlimits=sort([xlimits(1) xlimits(end)]);
ylimits=get(hi,'ydata');
ylimits=sort([ylimits(1) ylimits(end)]);

limdata=get(hlimbox,'userdata');
limlnpositions=limdata{2};
% [top bottom left side right side]
lntop=limlnpositions(1);
ytop=get(lntop,'ydata'); ytop=ytop(2);
lnbot=limlnpositions(2);
ybot=get(lnbot,'ydata'); ybot=ybot(1);
lnlft=limlnpositions(3);
xlft=get(lnlft,'xdata'); xlft=xlft(1);
lnrgt=limlnpositions(4);
xrgt=get(lnrgt,'xdata'); xrgt=xrgt(2);
[ytop ybot xlft xrgt];
limptpositions=limdata{1};
% this is cheating, but I don't want to move the points here, so I will
% just make their visibility 'off'... hehehe... tricky me
set(limptpositions,'visible','off');
% [upper left upper right lower left lower right]
limcent=limdata{4};
 
opt=get(axes1,'userdata');
xdat=get(axes1,'xlim');
ydat=get(axes1,'ylim');
cpt=get(axes1,'currentpoint');
newxpt=sort([xdat(1) cpt(1,1) xdat(2)]);
newypt=sort([ydat(1) cpt(1,2) ydat(2)]);
newpt=[newxpt(2) newypt(2)];    % limits for instant axes limits

% this will not let the cent out of the axis
set(limcent,'xdata',newpt(1),'ydata',newpt(2));

dx=opt(1,1)-cpt(1,1);
dy=opt(1,2)-cpt(1,2);
set(axes1,'userdata',cpt);
newxdat=sort([xlft-dx xrgt-dx]);
newydat=sort([ybot-dy ytop-dy]);

if(newxdat(1)<=xlimits(1)|newxdat(2)>=xlimits(2))
    newxdat=[xlft xrgt];
end
if(newydat(1)<=ylimits(1)|newydat(2)>=ylimits(2))
    newydat=[ybot ytop];
end
mdat=[newydat newxdat];
set(lntop,'xdata',[mdat(3) mdat(4)],'ydata',[mdat(1) mdat(1)]);
set(lnbot,'xdata',[mdat(3) mdat(4)],'ydata',[mdat(2) mdat(2)]);
set(lnlft,'xdata',[mdat(3) mdat(3)],'ydata',[mdat(2) mdat(1)]);
set(lnrgt,'xdata',[mdat(4) mdat(4)],'ydata',[mdat(2) mdat(1)]);


