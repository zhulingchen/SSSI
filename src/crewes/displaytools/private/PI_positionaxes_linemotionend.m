function PI_positionaxes_linemotionend(hObject, eventdata, handles)
global OLDDOWNFCN OLDMOTIONFCN OLDUPFCN
h=get(gcf,'userdata');
hmsg=h(2);
hi=h(5);
set(gcf,'windowbuttonmotionfcn',[]);
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
set(get(mainax,'children'),'erasemode','normal');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
lns=get(get(posax,'title'),'userdata');
set(lns(1),'erasemode','normal');
set(lns(2),'erasemode','normal');
set(lns(3),'erasemode','normal');
set(lns(4),'erasemode','normal');
set(lns(6),'erasemode','normal','visible','on');
h=get(gcf,'userdata');
hzoompick=h(9);
%restore button functions
set(gcf,'windowbuttondownfcn',OLDDOWNFCN);
set(gcf,'windowbuttonupfcn',OLDUPFCN);
set(gcf,'windowbuttonmotionfcn',OLDMOTIONFCN);
%value=get(hzoompick,'value');
% switch value
%     case 1
%         selboxinit('plotimage(''zoom'')',1);
%         set(gcf,'name','Seismic Image Plot, Simplezooming installed (Use MB1)');
%     case 2
%         drawlineinit('plotimage(''pick'')');
%         set(gcf,'name','Seismic Image Plot, Picking resummed (Use MB1)');
%     case 3
%         drawlineinit('plotimage(''pick'')');
%         set(gcf,'name','Seismic Image Plot, Picking new (Use MB1)');
%         set(hzoompick,'userdata',[]);
%     case 4
% %         udat=get(hi,'userdata');
% %         seisstruct.SEIS=udat{1};
% %         seisstruct.X=udat{2};
% %         seisstruct.T=udat{3};
% %         pref.hmsg=hmsg;
% %         masteraxes=gca;
% %         holdhandle=gca;
% %         picksle('plotimage(''pick'')',seisstruct,masteraxes,pref);
% %         set(gcf,'name','Seismic Image Plot, Picking resummed (Use MB1)');
% %         stringinfo='Automatic Picking has been enabled';
%         
%     
% end
ylim=sort([get(lns(1),'ydata') get(lns(2),'ydata')]);
xlim=sort([get(lns(3),'xdata') get(lns(4),'xdata')]);
set(mainax,'xlim',[xlim(1) xlim(end)],'ylim',[ylim(1) ylim(end)]);
set(hi,'alphadata',[1]);
PI_zoomlock;
