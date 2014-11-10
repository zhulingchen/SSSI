function PI_zoomlock(hObject, eventdata, handles)
global ZOOM_LOCKS
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
hgcf=gcf;
if(~isempty(ZOOM_LOCKS))
    hgca=mainax;
    ydat=get(hgca,'ylim');
    xdat=get(hgca,'xlim');
    for ii=1:size(ZOOM_LOCKS,1)
        CheckLock=ZOOM_LOCKS(ii,:);
        if(CheckLock(2)==hgcf)
            SlaveAxes=findobj(CheckLock(1),'type','axes','tag','MAINAXES');
            SlaveAxesChildren=get(SlaveAxes,'children');
            set(SlaveAxes,'ylim',ydat,'xlim',xdat);
            hslave=get(SlaveAxes,'Parent');
            set(0,'currentfigure',hslave);
            hslavedat=get(hslave,'userdata');
            hvertscrol=hslavedat(16);
            hhorscrol=hslavedat(17);
            % sliders on other figures are being shut off to avoid
            % problems with sliders being set wrong
            set(hhorscrol,'visible','off');
            set(hvertscrol,'visible','off');
            PI_positionaxes_lineposition;
        end
    end
end
set(0,'currentfigure',hgcf);
