function hbutt=raymigbutton(figseis,figmod,vel,dx,msg,clr)

if(isnumeric(figseis))
    action='init';
else
    action=figseis;
end
if(strcmp(action,'init'))
    if(nargin<5)
        msg='';
    end
    if(nargin<6)
        clr='r';
    end
    %see if there is already a button
    hkids=allchild(figseis);
    startpos=.02;
    for k=1:length(hkids)
        tag=get(hkids(k),'tag');
        if(strcmp(tag,'raymigbutton')||strcmp(tag,'clearraysbutton')||...
                strcmp(tag,'clearpicksbutton')||strcmp(tag,'raymodbutton'))
            pos=get(hkids(k),'position');
            rightedge=pos(1)+pos(3);
            if(rightedge>startpos)
                startpos=rightedge+.01;
            end
        end
    end
    if(isempty(msg))
        label=['Migrate fig' int2str(figseis) ' -> fig' int2str(figmod)];
        wpixels=120;%desired width in pixels
        %determine window size
        pos=get(gcf,'position');%should be in pixels
        width=wpixels/pos(3);%width in normaliized units
    else
        label=['Migrate fig' int2str(figseis) ' -> fig' int2str(figmod) ' ' msg];
        wpixels=150;%desired width in pixels
        %determine window size
        pos=get(gcf,'position');%should be in pixels
        width=wpixels/pos(3);%width in normaliized units
    end
    hbutt=uicontrol(figseis,'style','pushbutton','string',label,'units','normalized',...
        'position',[startpos .92 width .05],'callback','raymigbutton(''migrate'')',...
        'userdata',{vel;dx;figseis;figmod;clr},'tag','raymigbutton');
elseif(strcmp(action,'migrate'))
    ud=get(gcbo,'userdata');
    rayvelmod(ud{1},ud{2})
    eventraymig(ud{3},ud{4},nan,ud{5});
end
    