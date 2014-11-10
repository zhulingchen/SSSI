function PI_KillLimitLines(arg1)
% Killing limit lines
global SCALE_OPT NUMBER_OF_COLORS GRAY_PCT CLIP COLOR_MAP NOBRIGHTEN PICKS PICKCOLOR XAXISTOP ZOOM_VALUE ZOOM_LOCKS

if(nargin<=0)
else
    hfig=arg1;
    h=get(hfig,'userdata');
    hi=h(5);
    hlimbox=h(14);
    set(hlimbox,'userdata',[],'value',[1]);
    set(findobj(hfig,'type','uimenu','tag','LIMITBOXMENU'),'label','Limit Box: ON');
    for ii=1:size(PICKS,1)
        if(PICKS{ii,1}==hfig)
            PICKS{ii,2}=[];
            PICKS{ii,3}=[];
            break
        end
    end
    mfig=findobj(0,'type','figure','tag','MVLINESMEASUREMENTS');

    for ii=1:length(mfig)
        checkfig=get(mfig(ii),'userdata');
        if(checkfig==hfig)
            delete(checkfig);
            break
        end
    end
    delete(findobj(hfig,'tag','LIMITLINE'));
    delete(findobj(hfig,'tag','LIMITPOINT'));
    delete(findobj(hfig,'tag','PICKTEXT'));
    delete(findobj(hfig,'tag','PICKMARKER'));
    delete(get(hi,'uicontextmenu'));
end
