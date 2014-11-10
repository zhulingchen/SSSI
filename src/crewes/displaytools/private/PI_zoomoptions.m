function PI_zoomoptions;
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS

xdat=get(gca,'xlim');
ydat=get(gca,'ylim');
val=get(gcbo,'userdata');
switch val(2);
case 1
    % Publishing current axes zoom values
    ZOOM_VALUE=[xdat(1) xdat(2) ydat(1) ydat(2) gcf];
case 2
    % Setting current axes to published zoom values
    if(~isempty(ZOOM_VALUE))
        set(gca,'xlim',sort([ZOOM_VALUE(1:2)]),'ylim',sort([ZOOM_VALUE(3:4)]));
    end
case 3
    % locking current plot to another plot
    FigureNumber=val(1);
    FigureAxes=findobj(FigureNumber,'type','axes','tag','MAINAXES');
    xdat=get(FigureAxes,'xlim');
    ydat=get(FigureAxes,'ylim');
    set(gca,'xlim',xdat,'ylim',ydat);
    for ii=1:size(ZOOM_LOCKS,1);
        xx=ZOOM_LOCKS(ii,:);
        CheckLock=find(xx==gca);
        if(~isempty(CheckLock))
        end
    end
    ZOOM_LOCKS=[ZOOM_LOCKS;gcf FigureNumber];
case 4
    % unlocking current plot
    for ii=1:size(ZOOM_LOCKS,2)
        CheckLock=PI_ZOOM_LOCKS(ii);
        if(CheckLock(1)==gcf)
            ZOOM_LOCKS(ii,:)=[];
        end
    end
end
