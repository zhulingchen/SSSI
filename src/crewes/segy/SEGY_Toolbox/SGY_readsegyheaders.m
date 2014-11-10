function [texthead,binaryhead,tracehead,extendedhead]=SGY_readsegyheaders(filein)
texthead=TextHeader(filein);
binaryhead=BinaryHeader(filein);
extendedhead={};
%test to see if there is an extended header
fseek(binaryhead.fid,3504,'bof');
numofext=fread(binaryhead.fid,1,'uint16');
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
            extendedhead{n,1}=TextHeader(filein,'txthoffset',num2str(trstat));
        else
            break;
        end
        trstat=trstat+3200;
    end
    
end

tracehead=TraceHeader(filein,'hdroffset',num2str(trstat),'extendedheaderflag','no');
end
