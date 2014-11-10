function hbutt=raymodbutton(figmod,figseis,vel,dx,msg,clr)

if(isnumeric(figmod))
    action='init';
else
    action=figmod;
end
if(strcmp(action,'init'))
    if(nargin<5)
        msg='';
    end
    if(nargin<6)
        clr='r';
    end
    %see if there is already a button
    hkids=allchild(figmod);
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
        label=['Model fig' int2str(figmod) ' -> fig' int2str(figseis)];
        wpixels=120;%desired width in pixels
        %determine window size
        pos=get(gcf,'position');%should be in pixels
        width=wpixels/pos(3);%width in normaliized units
    else
        label=['Model fig' int2str(figmod) ' -> fig' int2str(figseis) ' ' msg];
        wpixels=150;%desired width in pixels
        %determine window size
        pos=get(gcf,'position');%should be in pixels
        width=wpixels/pos(3);%width in normaliized units
    end
    hbutt=uicontrol(figmod,'style','pushbutton','string',label,'units','normalized',...
        'position',[startpos .92 width .05],'callback','raymodbutton(''model'')',...
        'userdata',{vel;dx;figmod;figseis;clr},'tag','raymodbutton');
elseif(strcmp(action,'model'))
    ud=get(gcbo,'userdata');
    rayvelmod(ud{1},ud{2})
    eventraymod(ud{3},ud{4},nan,ud{5});
end
    