function PI_lmptreset(action);
h=get(gcf,'userdata');
if(strcmp(action,'lmptreset'))
    slavefig=gcf;
    hmaster=h(6);
    set(0,'currentfigure',hmaster');
    h=get(gcf,'userdata');
else
    hlimbox=h(14);
    limdat=get(hlimbox,'userdata');
    slavefig=limdat{3};
end
hmsg=h(2);
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
limlndat=limdat{2};
% [top bottom left right] - reseting Lines
parentaxis=get(limlndat(1),'parent');
xdat=get(parentaxis,'xlim');
ydat=get(parentaxis,'ylim');
set(limlndat(1),'ydata',[ydat(2) ydat(2)],'xdata',[xdat(1) xdat(2)]);
set(limlndat(2),'ydata',[ydat(1) ydat(1)],'xdata',[xdat(1) xdat(2)]);
set(limlndat(3),'xdata',[xdat(1) xdat(1)],'ydata',[ydat(1) ydat(2)]);
set(limlndat(4),'xdata',[xdat(2) xdat(2)],'ydata',[ydat(1) ydat(2)]);
allnums=[ydat xdat];
findfig=limdat{3};
mlninfo=get(findfig,'userdata');
ddat={'ydata' 'ydata' 'xdata' 'xdata'};
for ii=1:4
    set(mlninfo(ii),'string',num2str(allnums(ii)));
end
limptdat=limdat{1};
% [upperleft upperright lowerleft lowerright] - reseting points
set(limptdat(1),'ydata',ydat(2),'xdata',xdat(1));
set(limptdat(2),'ydata',ydat(2),'xdata',xdat(2));
set(limptdat(3),'ydata',ydat(1),'xdata',xdat(1));
set(limptdat(4),'ydata',ydat(1),'xdata',xdat(2));
plotimage('limboxrescale');

limcent=limdat{4};
set(limcent,'ydata',(ydat(2)-ydat(1))/2+ydat(1),'xdata',(xdat(2)-xdat(1))/2+xdat(1));
stringinfo='Limit lines have been reset to match axes.';
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
