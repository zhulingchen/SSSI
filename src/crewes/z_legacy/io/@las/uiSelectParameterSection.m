function varargout = uiSelectParameterSection(obj, section, msg)
%
% function varargout = uiselectsection(obj, section, msg)
%
% UISELECTSECTION is a GUI-based utility designed to allow the selection of
% sections in LAS files. This is useful for files that contain multiple
% data sections
%

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its author (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

if nargin < 2 || nargin > 3
    warning('crewes:las:uiselectparamatersection','Not enough input parameters');
    return
elseif isequal(nargin,2)
    msg = '';
end

%assume failure
varargout{1} = '';
varargout{2} = [];

hfig=figure('name','las.uiselectparametersection','numbertitle','off','units','normalized','position',[.3 .3 .6 .4]);
set(hfig,'menubar','none');
set(hfig,'closerequestfcn',@closereq_callback);
    
if iscell(section)
    n = numel(section{1});
    if ~isequal(n,4)
        warning('crewes:las:uiselectsection',...
            'Input cell array ''section'' does appear to be las object format');
        return
    else
        sidx = 1:length(section);
        dts = cell2mat(cellfun(@(x) x{4},section,'UniformOutput',false));
        ldat = cellfun(@(x) x{1},section,'UniformOutput',false)';
        ldat = [ldat{:}];
        lmnem = cellfun(@(x) x{3}(1,:),section,'UniformOutput',false)';
   end
elseif ischar(section)
    sidx = obj.getSectionIndices(section);
    dts  = obj.isDataSection(sidx);
    ldat  = obj.getSectionNames(section)';
    lmnem = obj.getParameterSectionMnemonics(ldat)';
end

if numel(sidx) < 2
    delete(hfig)
    if ~dts
        varargout{1} = lmnem;
        varargout{2} = sidx;
    end
    return
end

if isempty(sidx) || sum(~dts) ~= numel(dts)
    delete(hfig)
    warning('crewes:las:uiselectsection',...
        'Section name does not exclusively match parameter sections');
    return
end

tabledata = cell(numel(ldat),3);
tabledata(:,1)={false};
tabledata(1,1)={true};

for n = 1:numel(ldat)
    tabledata(n,2) = ldat(n);
    tabledata(n,3) = {sprintf('%s ',lmnem{n}{:})};
end

deswidth = max(cellfun('length',tabledata(:,3)));
daswidth = max(cellfun('length',tabledata(:,2)));
das = 'Parameter Section';
des = 'Mnemonics in section';
das = [das char(32*ones(1,2*daswidth-length(das)))];
des = [des char(32*ones(1,2*deswidth-length(des)))];

colnames={'Select',das,des};
colformat={'logical','char','char'};
columneditable =  [true false false];
posn=[.1 .1 .7 .8];
uitable('Parent',hfig,'Tag','logs_table','data',tabledata,'columnname',colnames,...
    'columnformat',colformat,'columneditable',columneditable,...
    'rowname',[],'units','normalized','position',posn,...
    'celleditcallback',@celledit_callback);

%add some buttons
xnow=.8;ynow=.8;width=.1;height=.1;
uicontrol('Parent',hfig,'Tag','next_button','style','pushbutton','string','Next...',...
    'callback',@next_callback,'visible','on',...
    'units','normalized','position',[xnow,ynow,width,height]);

xnow=.8;ynow=.7;width=.1;height=.1;
uicontrol('Parent',hfig,'Tag','cancel_button','style','pushbutton','string','Cancel',...
    'callback',@closereq_callback,'visible','on',...
    'units','normalized','position',[xnow,ynow,width,height]);

xnow=.8;ynow=.1;width=.1;height=.6;
uicontrol('Parent',hfig,'Tag','msg_text','style','text','string',msg,...
    'visible','on','units','normalized','position',[xnow,ynow,width,height])

%put a title at the top
xnow=.1;ynow=.95;width=.7;height=.05;
uicontrol('Parent',hfig,'Tag','filename_text','style','text','string',['File:' ' ' obj.fileName],...
    'visible','on',...
    'units','normalized','position',[xnow,ynow,width,height]);
xnow=.1;ynow=.9;width=.7;height=.05;

if str2double(obj.version) < 3.0
    well_parameter = '~w';
else
    well_parameter = '~well';
end
name   = obj.getMnemonicValue(well_parameter,'WELL');
wellid = obj.getMnemonicValue(well_parameter,'UWI');
uicontrol('Parent',hfig,'Tag','wellid_text','style','text','string',[name ' ' wellid],...
    'visible','on',...
    'units','normalized','position',[xnow,ynow,width,height]);
xnow=.1;ynow=.05;width=.7;height=.05;

uicontrol('Parent',hfig,'Tag','selected_section','style','text','string',...
    ['Section: ' ldat{1} ' selected for editing'],...
    'visible','on',...
    'units','normalized','position',[xnow,ynow,width,height]);

uicontrol('Parent',hfig,'Tag','done_text','style','text','string',...
    '',...
    'visible','off',...
    'units','normalized','position',[xnow,ynow,width,height]);

handles = guihandles(hfig); %create struct based on 'Tag' property
guidata(hfig,handles);      %save struct in hfig

set(handles.done_text,'UserData',0);
set(handles.selected_section,'UserData', {ldat{1},1});

waitfor(handles.done_text,'UserData');

t = get(handles.done_text,'UserData');
switch t
    case 1 %cancel button clicked, or figure close
        warning('crewes:las:uiselectsection',...
             ['User canceled section selection; '...
              'Keeping first section in LAS file']);
        varargout = get(handles.selected_section,'UserData');
    case 2 %section selected
        varargout = get(handles.selected_section,'UserData');
end
delete(hfig);

end %end function

function closereq_callback(varargin)
handles = guidata(gcbo);
set(handles.done_text,'UserData',1);
guidata(gcbo);
end %end function closereq_callback

function celledit_callback(varargin)

handles = guidata(gcbo);
d = get(handles.logs_table,'data');

events = varargin{2};

d(:,events.Indices(2))={false}; %Set all cells in checkbox column to false
d(events.Indices(1),events.Indices(2))={events.EditData}; %Set cell clicked on to editdata
set(handles.logs_table,'data',d);

switch (events.EditData)
    case 1,
        set(handles.selected_section,'string',...
            ['Section: ', d{events.Indices(1),2}, ...
            ' selected for editing']);
        set(handles.selected_section,'UserData', {d{events.Indices(1),2}...
            events.Indices(1)} );
    case 0,
        set(handles.selected_section,'string',...
            'No sections selected for editing');
        set(handles.selected_section,'UserData', cell(1,2));
    otherwise,
        set(handles.selected_section,'UserData', cell(1,2));
end

guidata(gcbo);
end %end function celledit_callback

function next_callback(varargin)
handles = guidata(gcbo);

section2edit = get(handles.selected_section,'UserData');
set(handles.next_button,'UserData',section2edit);
set(handles.done_text,'UserData',2);
guidata(gcbo);

end %end function next_callback
