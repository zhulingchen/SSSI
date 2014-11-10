function traceheaderviewer(tracehead)

try
if ~isa(tracehead,'TraceHeader')
    me=MExecption('traceheaderviewer:InvalidInputType',...
        'tracehead must be a TraceHeader Object');
    throw(me)
end

indshnm=strcmp(tracehead.definitions.keys,'Name');
inddesc=strcmp(tracehead.definitions.keys,'Description');

dat={};
for k=1:length(tracehead.definitions.values(:,indshnm))
    rowname{k,1}=tracehead.definitions.values{k,indshnm};
    dat{k,1}=tracehead.definitions.values{k,inddesc};
    dat{k,2}=num2str(unique(getfield(tracehead.header,tracehead.definitions.values{k,indshnm})));
end

h.fig=figure('menubar','none','numbertitle','off','name','Trace Header Viewer',...
    'units','normalized','position',[.1,.25,.8,.5]);
h.table=uitable(h.fig,'units','normalize','position',[.05,.15,.9,.8],...
    'columnname',{'Description','Possible Values'},...
    'columnwidth',{400, 500},'RowName',rowname,'data',dat,...
    'CellSelectionCallback',@seecellcontents );
h.ok=uicontrol('style','pushbutton','string','OK','units','normalized',...
    'Position',[.45,.04,.1,.07],'callback',@selectok);
guidata(h.fig,h);



catch me
    error(me.message,me.identifier);
end

    function seecellcontents(hObject,eventdata)
        h=guidata(hObject);
        dat=get(h.table,'data');
        str=dat{eventdata.Indices(1),eventdata.Indices(2)};
        msgbox(str);
    end

    function selectok(hObject, eventdata)
        h=guidata(hObject);
        delete(h.fig);
    end
        



end