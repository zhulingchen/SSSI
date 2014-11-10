function helpbutton(figno,helpmsg,titlstr)
% Add a help button to a plotimage figure
if(isnumeric(figno))
    action='init';
else
    action=figno;
end

if(strcmp(action,'init'))
    %check for existing buttons
    hkids=allchild(figno);
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
    wpixels=60;%desired width in pixels
    %determine window size
    pos=get(gcf,'position');%should be in pixels
    %
    width=wpixels/pos(3);%width in normaliized units
    uicontrol(figno,'style','pushbutton','string','Help!','units','normalized',...
        'position',[startpos .92 width .05],'callback','helpbutton(''help'')',...
        'tag','clearpicksbutton','userdata',{helpmsg;titlstr});
elseif(strcmp(action,'help'))
    hbutt=gcbo;
    udat=get(hbutt,'userdata');
    helpmsg=udat{1};
    titlestring=udat{2};
    msgbox(helpmsg,titlestring);
end