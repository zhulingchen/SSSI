function segyobj=readsegyheaders(segyobj,filein,varargin)
segyobj.textheader=TextHeader(filein,varargin{:});
segyobj.binaryheader=BinaryHeader(filein,varargin{:});
segyobj.extendedheader={};
%test to see if there is an extended header
if ~strcmpi(varargin,'extendedheader')
fseek(segyobj.fid,3504,'bof');
numofext=fread(segyobj.fid,1,'uint16');
trstat=3600;
screensz=get(0,'ScreenSize');
if numofext
    if numofext==-1
        numofext=500;
    end
    for n=1:numofext
        extheadtmp=TextHeader(filein,'txthoffset',num2str(trstat));
        response=listdlg('PromptString','Does this appear to be a somewhat legible Extended Header?',...
            'ListString',extheadtmp.header,'Name','Is This An Extended Header?',...
            'OkString','Yes','CancelString','No','SelectionMode','Single',...
            'ListSize',[.5*screensz(3),.5*screensz(4)]);
        if isempty(response)
            break
        else
            if response
                segyobj.extendedheader{n,1}=TextHeader(filein,'txthoffset',num2str(trstat));
                trstat=trstat+3200;
            else
                break;
            end
        end
    end
    
end

segyobj.traceheader=TraceHeader(filein,'hdroffset',num2str(trstat),'extendedheaderflag','no');
else
    if strcmpi(varargin,'hdroffset');
        ind=strfind(varargin,'hdroffset');
        segyobj.traceheader=TraceHeader(filein,'hdroffset',varargin{ind},'extendedheaderflag','no');
    else
        segyobj.traceheader=TraceHeader(filein,'extendedheaderflag','no');
    end
end
end
