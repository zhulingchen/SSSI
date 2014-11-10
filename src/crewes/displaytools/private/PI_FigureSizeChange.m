function PI_FigureSizeChange
% changes the look of the plotimage figure for publishable purposes
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
haxs=findobj(gcf,'type','axes','tag','MAINAXES');
h=get(gcf,'userdata');
hmsg=h(2);
ttl=get(haxs,'title');
dat=get(ttl,'userdata');
if(isempty(dat))
    return
end
checkcolour=get(gcf,'color');
if(checkcolour~=[1 1 1]);
    dat{4}=get(gcf,'position');
    set(ttl,'userdata',dat);
    whitefig;
    bigfont(gcf,1.5);
    boldlines(gcf,2,2);
    bigfig;
    stringinfo=['Publishable Plotimage'];
    set(hmsg,'visible','off');
    set(gcbo,'label','Publishable Figure: ON');
else
    siz=dat{4};
    greyfig;
    bigfont(gcf,1/1.5,1);
    boldlines(gcf,.5,.5);
    set(gcf,'position',siz);
    stringinfo=['Normal Plotimage'];
    set(hmsg,'visible','on');
    set(gcbo,'label','Publishable Figure: OFF');
end
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
