function PI_limlnmove2();
hobj=gco;
axes1=gca;
pt=get(axes1,'userdata');
xdat=get(axes1,'xlim');
ydat=get(axes1,'ylim');
cpt=get(axes1,'currentpoint');
newxln=sort([xdat(1) cpt(1,1) xdat(2)]);
newyln=sort([ydat(1) cpt(1,2) ydat(2)]);
newlndata=[newyln(2) newxln(2)];
% getting all line data to find out which line is being moved
% and move others depending on location
h=get(gcf,'userdata');
hlimbox=h(14);
limdata=get(hlimbox,'userdata');
limcent=limdata{4};
limlnpositions=limdata{2};
lntop=limlnpositions(1);
lnlft=limlnpositions(3);
limptpositions=limdata{1};
% [top bottom left side right side] - only original configuration
checkln=(find(limlnpositions==hobj));
shiftcell1={'ydata' 'ydata' 'xdata' 'xdata'};
% Shiftcell1 Configuration
% 1 / 2 - Lines perpendicular to Line being moved
% 3 - Line parrallel to line being moved
% 4 / 5 - Points on line that are being moved
shiftcell2={[3 4 2 1 2] [3 4 1 3 4] [1 2 4 1 3] [1 2 3 2 4]};
shiftpos=shiftcell2{checkln};
set(hobj,shiftcell1{checkln},[newlndata(round(checkln/2)) newlndata(round(checkln/2))]);
for ii=1:2
    xlnx=shiftpos(ii);
    shiftln=limlnpositions(xlnx);
    xshiftx=shiftcell2{checkln};
    holdln=limlnpositions(xshiftx(3));
    holdata=get(holdln,shiftcell1{checkln});
    newposdata=sort([holdata(1) newlndata(round(checkln/2))]);
    set(shiftln,shiftcell1{checkln},[newposdata(1) newposdata(2)]);
    % moving first points
    shiftpt=limptpositions(xshiftx(ii+3));
    set(shiftpt,shiftcell1{checkln},newlndata(round(checkln/2)));
end
xdat=get(lntop,'xdata');
ydat=get(lnlft,'ydata');
set(limcent,'xdata',(xdat(2)-xdat(1))/2+xdat(1),'ydata',(ydat(2)-ydat(1))/2+ydat(1));
