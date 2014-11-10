function PI_MovePickLineEnd();
global PICKS 
h=get(gcf,'userdata');
hmsg=h(2);
set(findobj(gcf,'type','line','tag','PICKS'),'erasemode','normal');
set(gcf,'windowbuttonmotionfcn','','windowbuttonupfcn','');
hobj=gco;
nm=get(hobj,'tag');
udat=get(hobj,'userdata');
switch nm
    case 'PICKMARKER'
        set(hobj,'erasemode','normal');
        set(udat{1},'erasemode','normal');
        xdat=get(udat{1},'xdata');
        ydat=get(udat{1},'ydata');
        hpick=udat{1};
    case 'PICKS'
        xdat=get(hobj,'xdata');
        ydat=get(hobj,'ydata');
        set(udat(1),'erasemode','normal','ydata',ydat(1),'xdata',xdat(1),...
            'visible','on');
        set(udat(2),'erasemode','normal','ydata',ydat(2),'xdata',xdat(2),...
            'visible','on');
        hpick=hobj;
end
if(~isempty(PICKS))
    for ii=1:size(PICKS,1)
        CheckFigure=PICKS{ii,1};
        if(CheckFigure==gcf)
            CheckHandles=PICKS{ii,3};
            for jj=1:length(CheckHandles)
                if(CheckHandles(jj)==hpick)
                    LinePositions=PICKS{ii,2};
                    LinePositions(jj,:)=[xdat(1) ydat(1) xdat(2) ydat(2)];
                    PICKS{ii,2}=LinePositions;
                end
            end
        end
    end
end
set(findobj(gcf,'type','line','tag','LIMITLINE'),'erasemode','normal');
set(findobj(gcf,'type','line','tag','LIMITPOINT'),'erasemode','normal');
set(findobj(gcf,'type','text','tag','PICKTEXT'),'erasemode','normal','visible','off');
set(gcf,'windowbuttonmotionfcn','','windowbuttonupfcn','');
stringinfo='Pick line has been moved';
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
