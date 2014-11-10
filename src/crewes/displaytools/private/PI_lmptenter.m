function PI_lmptenter();
hmeasurefig=get(gcf,'userdata');
hmaster=hmeasurefig(6);
set(0,'currentfigure',hmaster');
h=get(gcf,'userdata');
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
lmptdat=limdat{1};
lmlndat=limdat{2};
lmcent=limdat{4};
parentaxis=get(lmptdat(1),'parent');
% ensuring that axis does not change dimensions
xdat=get(parentaxis,'xlim');
ydat=get(parentaxis,'ylim');
ckmdat=[ydat xdat];
mdat=[];
% [Bottom top left right]
for ii=1:4
    nm=str2num(get(hmeasurefig(ii),'string'));
    mdat=[mdat nm];
end
ckydat=sort([ckmdat(1) ckmdat(2) mdat(1) mdat(2)]);
ckxdat=sort([ckmdat(3) ckmdat(4) mdat(3) mdat(4)]);
mdat=[ckydat(2) ckydat(3) ckxdat(2) ckxdat(3)];
% [top bottom left right] - reseting Lines
hold on;
set(lmlndat(1),'ydata',[mdat(2) mdat(2)],'xdata',[mdat(3) mdat(4)]);
set(lmlndat(2),'ydata',[mdat(1) mdat(1)],'xdata',[mdat(3) mdat(4)]);
set(lmlndat(3),'ydata',[mdat(1) mdat(2)],'xdata',[mdat(3) mdat(3)]);
set(lmlndat(4),'ydata',[mdat(1) mdat(2)],'xdata',[mdat(4) mdat(4)]);
% [upperleft upperright lowerleft lowerright] - reseting points
set(lmptdat(1),'ydata',mdat(2),'xdata',mdat(3));
set(lmptdat(2),'ydata',mdat(2),'xdata',mdat(4));
set(lmptdat(3),'ydata',mdat(1),'xdata',mdat(3));
set(lmptdat(4),'ydata',mdat(1),'xdata',mdat(4));
% setting measurement  values
for ii=1:4
    set(hmeasurefig(ii),'string',num2str(mdat(ii)));
end
set(lmcent,'ydata',(mdat(2)-mdat(1))/2+mdat(1),'xdata',(mdat(4)-mdat(3))/2+mdat(3));
set(parentaxis,'ylim',ydat,'xlim',xdat);
hold off;
plotimage('limboxrescale');
