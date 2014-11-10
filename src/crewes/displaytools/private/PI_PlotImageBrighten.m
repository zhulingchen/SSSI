function PI_PlotImageBrighten();
h=get(gcf,'userdata');
hlabel=h(11);
hslider=h(12);
currlvl=get(hlabel,'userdata');
newlvl=get(hslider,'value');
if(currlvl>newlvl)
    nsteps=currlvl-newlvl;
    for k=1:nsteps
        brighten(-.25);
    end
elseif(currlvl<newlvl)
    nsteps=newlvl-currlvl;
    for k=1:nsteps
        brighten(.25);
    end
end
set(hlabel,'string',['Bright ' int2str(newlvl)])
set(hlabel,'userdata',newlvl);
plotimage('limboxrescale');
