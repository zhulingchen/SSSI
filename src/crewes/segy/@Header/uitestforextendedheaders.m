function traceheaderposition=uitestforextendedheaders(obj)
fid=obj.fid;
filein=obj.filename;
trstat=obj.hdroffset;
traceheaderposition=obj.hdroffset;
if nargin<2
    trstat=3600;
end



fseek(fid,3504,'bof');
numofext=fread(fid,1,'uint16');
if numofext
    if numofext==-1
        numofext=500;
    end
    screensz=get(0,'ScreenSize');
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
                trstat=trstat+3200;
            else
                break;
            end
        end
        
    end
    traceheaderposition=trstat;
end
