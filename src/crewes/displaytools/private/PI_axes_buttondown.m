function PI_axes_buttondown
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
h=get(gcf,'userdata');
hzoompick=h(9);
value=get(hzoompick,'value');
if(mainax==gca)
    switch value
        case 1
            selboxinit('plotimage(''zoom'')',1);
        case 2
            selboxinit('plotimage(''zoom'')',1)
        case 3
            
    end
    
elseif(posax==gca)
    
end
