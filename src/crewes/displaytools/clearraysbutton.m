function clearraysbutton(figno)
% Add a clear rays button to a plotimage fighure

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
uicontrol(figno,'style','pushbutton','string','clear rays','units','normalized',...
    'position',[startpos .92 width .05],'callback','clearrays','tag','clearraysbutton');