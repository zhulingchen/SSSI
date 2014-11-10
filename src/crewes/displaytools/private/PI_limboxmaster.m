function PI_limboxmaster();
h=get(gcf,'userdata');
hi=h(5);
hmaster=h(10);
hlimbox=h(14);
limfig=get(hlimbox,'userdata');
limdat=get(limfig{3},'userdata');
limmat=[];
for ii=1:4
    limmat=[limmat str2num(get(limdat(ii),'string'))];
end
ydat=sort([round(limmat(1)) round(limmat(2))]);
xdat=sort([round(limmat(3)) round(limmat(4))]);
% checking to see if lines numbers are acceptable, right now values
% will not change if the new lines are within 3 % of eachother
if((abs(ydat(1)-ydat(2)))/(ydat(1)+ydat(2))*100<=3)|((abs(xdat(1)-xdat(2)))/(xdat(1)+xdat(2))*100<=3)
    return
end
seis=get(hi,'cdata');
seis2=seis(ydat(1):ydat(2),xdat(1):xdat(2));
mxs=full(max(max(abs(seis2))));
smean=full(mean(mean(seis2)));
stddev=full(sqrt( sum(sum((seis-smean).^2 ) )...
    /prod(size(seis))));
% set(hmaster,'userdata',[mxs smean stddev]);
% plotimage('rescale');
