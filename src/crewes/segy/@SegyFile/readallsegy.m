function segyobj=readallsegy(segyobj,filein,varargin)
% SGY_readallsegy(filein)
%
% SGY_readallsegy reads in the entire SEG-Y File 
name={};
value={};
 for i = 1:2:length(varargin)
                name = varargin{i};
                value = varargin{i+1};
 end

segyobj.textheader=TextHeader(filein);
segyobj.binaryheader=BinaryHeader(filein);

%test to see if there is an extended header
fseek(segyobj.binaryheader.fid,3504,'bof');
numofext=fread(segyobj.binaryheader.fid,1,'uint16');
trstat=3600;
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
                segyobj.extendedheader{n,1}=TextHeader(filein,'txthoffset',num2str(trstat));
                trstat=trstat+3200;
            else
                break;
            end
        end
        
    end
    
end
trch=Trace(filein,'hdroffset',num2str(trstat),'extendedheaderflag','no','new','1');
if any(strcmpi(name,'traces'))
    segyobj.trace=trch;
else
segyobj.trace=trch.getTraces();   
end

end
