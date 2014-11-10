function PI_zoomscroll();
h=get(gcf,'userdata');
hi=h(5);
imxdat=get(hi,'xdata');
imydat=get(hi,'ydata');
ha=gca;
axdat=get(ha,'xlim');
aydat=get(ha,'ylim');
delete(findobj(gcf,'type','line','tag','PICKMARKER'));
set(findobj(gcf,'type','line','tag','PICKS'),'erasemode','xor');
set(gcf,'pointer','fleur')
set(gcf,'windowbuttondownfcn','plotimage(''zoominout'')',...
    'windowbuttonmotionfcn','plotimage(''zoomscrollmotion'')',...
    'windowbuttonupfcn','plotimage(''zoomfcnend'')');
axspts=get(gca,'currentpoint');
figpts=get(gcf,'currentpoint');
set(gca,'userdata',[axspts(1,1) axspts(1,2) figpts(2)]);
set(hi,'alphadata',[75]);
return
