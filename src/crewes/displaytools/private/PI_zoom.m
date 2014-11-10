function PI_zoom;
global ZOOM_LOCKS
if(strcmp(get(gca,'visible'),'on'))
    set(gcf,'pointer','arrow');
    h=get(gcf,'userdata');
    hi=h(5);
    haxs=gca;
    hvertscrol=h(16);
    hhorscrol=h(17);
    % making sure that user is trying to zoom and not clicking on
    % on the limit box lines of points
    box=selboxfini;
    hi=h(5);
    if(length(box)==5)
        if(~ishandle(box(5)))
            return
        end
        delete(box(5));
    end
    if(isempty(box)|iscell(box)|length(box)<=3) return; end
    xmin=min([box(1) box(3)]);
    xmax=max([box(1) box(3)]);
    ymin=min([box(2) box(4)]);
    ymax=max([box(2) box(4)]);
    %get the current axis settings
    xlim=get(gca,'xlim');
    ylim=get(gca,'ylim');
    test1=xmin-xlim(1)+xmax-xlim(2)+ymin-ylim(1)+ymax-ylim(2);
    test2=(xmin-xmax)*(ymin-ymax);
    if(abs(test1)<10*eps|abs(test2)< 10*eps|strcmp(get(gcf,'selectiontype'),'extend'))
        xdat=get(hi,'xdata');
        ydat=get(hi,'ydata');
        xdat=[min(xdat) max(xdat)];
        ydat=[min(ydat) max(ydat)];
        set(gca,'xlim',xdat,'ylim',ydat);
        set(hhorscrol,'visible','off');
        set(hvertscrol,'visible','off');
        posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
        lns=get(get(posax,'title'),'userdata');
        set(lns(1),'visible','off');
        set(lns(2),'visible','off');
        set(lns(3),'visible','off');
        set(lns(4),'visible','off');
        set(lns(6),'visible','off');
    else
        imxdat=get(hi,'xdata');
        imydat=get(hi,'ydata');
        xdat=sort([imxdat(1) xmin xmax imxdat(end)]);
        ydat=sort([imydat(1) ymin ymax imydat(end)]);
        xdat=[xdat(2) xdat(3)];
        ydat=[ydat(2) ydat(3)];
        set(gca,'xlim',[xdat],'ylim',[ydat]);
        ximlim=sort(get(hi,'xdata'));
        yimlim=sort(get(hi,'ydata'));
        xdat2=sort(get(gca,'xlim'));
        ydat2=sort(get(gca,'ylim'));
        set(hhorscrol,'value',(xdat2(2)-xdat2(1))/2+xdat2(1),'visible','on','enable','on',...
            'userdata',{[2] [(xdat2(2)-xdat2(1))/2+xdat2(1)] []},...
            'max',ximlim(end),'min',ximlim(1));
        y1=yimlim(end)-(ydat2(1)-yimlim(1));
        y2=yimlim(1)+(yimlim(end)-ydat2(2));
        set(hvertscrol,'value',y2+(y1-y2)/2,'visible','on','userdata',{[1] [y2+(y1-y2)/2] []},...
            'enable','on','max',yimlim(end),'min',yimlim(1));
        PI_positionaxes_lineposition;
    end
    PI_zoomlock;   
else
end
set(findobj(gcf,'type','line','tag','LIMITLINE'),'erasemode','normal');
set(findobj(gcf,'type','line','tag','LIMITPOINT'),'erasemode','normal');
