function PI_zoomfcnend();
global ZOOM_LOCKS
set(findobj(gcf,'type','line','tag','PICKS'),'erasemode','normal');
set(gcf,'pointer','arrow');
set(gcf,'windowbuttonupfcn',[]);
set(gcf,'windowbuttonmotionfcn',[]);
if(~isempty(ZOOM_LOCKS))
    for ii=1:size(ZOOM_LOCKS,1)
        CheckLock=ZOOM_LOCKS(ii,:);
        if(CheckLock(2)==gcf)
            SlaveAxes=get(CheckLock(1),'currentaxes');
            SlaveAxesChildren=get(SlaveAxes,'children');
            set(SlaveAxesChildren,'erasemode','normal');
        end
    end
end
selboxinit('plotimage(''zoom'')',1);
h=get(gcf,'userdata');
hi=h(5);
set(hi,'alphadata',[1]);
