function PI_limlnfigurevis(action);
if(strcmp(action,'limlnfigurevis2'))
    % visiblilty being shut off by measurement figure
    set(gcf,'visible','off');
    h=get(gcf,'userdata');
    hvisbutton=h(5);
else
    h=get(gcf,'userdata');
    hmsg=h(2);
    hlimbox=h(14);
    limdat=get(hlimbox,'userdata');
    slavefig=limdat{3};
    j=get(slavefig,'userdata');
    jvisbutton=j(5);
    val=get(jvisbutton,'value');
    vis={'ON' 'OFF'};
    vvl=[2 1];
    set(slavefig,'visible',vis{val});
    set(jvisbutton,'value',vvl(val));
    set(hmsg,'string',[ 'Setting Limit Line figure visibility: ' vis{val}],'backgroundcolor',[1 1 1]);
end
