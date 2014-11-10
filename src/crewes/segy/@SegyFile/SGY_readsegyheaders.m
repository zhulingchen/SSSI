function segyobj=SGY_readsegyheaders(segyobj,filein)
segyobj.textheader=TextHeader(filein);
segyobj.binaryheader=BinaryHeader(filein);
segyobj.extendedheader={};
%test to see if there is an extended header
fseek(segyobj.fid,3504,'bof');
numofext=fread(segyobj.fid,1,'uint16');
trstat=3600;
if numofext
    if numofext==-1
        numofext=500;
    end
    for n=1:numofext
        extheadtmp=TextHeader(filein,'txthoffset',num2str(trstat));
        response=questdlg(extheadtmp.header,'Is This An Extended Header?','yes','No','No');
        if isempty(response)
            break
        elseif strcmp(response,'yes')
            segyobj.extendedheader{n,1}=TextHeader(filein,'txthoffset',num2str(trstat));
        else
            break;
        end
        trstat=trstat+3200;
    end
    
end

segyobj.traceheader=TraceHeader(filein,'hdroffset',num2str(trstat),'extendedheaderflag','no');
end
