function PI_positionaxes_linebuttondown(hObject, eventdata, handles)
global OLDDOWNFCN OLDMOTIONFCN OLDUPFCN
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
set(get(mainax,'children'),'erasemode','xor');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
lns=get(get(posax,'title'),'userdata');
set(lns(1),'erasemode','xor');
set(lns(2),'erasemode','xor');
set(lns(3),'erasemode','xor');
set(lns(4),'erasemode','xor');
set(lns(6),'erasemode','xor');
h=get(gcf,'userdata');
hmsg=h(2);
h=get(get(posax,'title'),'userdata');
pt=get(gca,'currentpoint');
set(gca,'userdata',pt);
delete(findobj(gcf,'type','line','tag','PICKMARKER'));
delete(findobj(gcf,'type','text','tag','PICKTEXT'));
OLDDOWNFCN=get(gcf,'windowbuttondownfcn');
OLDMOTIONFCN=get(gcf,'windowbuttonmotionfcn');
OLDUPFCN=get(gcf,'windowbuttonupfcn');
set(gcf,'windowbuttonupfcn',@PI_positionaxes_linemotionend);
set(gcf,'windowbuttonmotionfcn',@PI_positionaxes_linemotion);
