function [texthead,binaryhead,traceh,traces,extendedhead]=SGY_readallsegy(filein)
% SGY_readallsegy(filein)
%
% SGY_readallsegy reads in the entire 


texthead=TextHeader(filein);
binaryhead=BinaryHeader(filein);

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
        if strcmp(response,'yes')
            extendedhead{n,1}=TextHeader(filein,'txthoffset',num2str(trstat));
        else
            break;
        end
        trstat=trstat+3200;
    end
    
end
trch=TraceHeader(filein,'hdroffset',num2str(trstat),'extendedheaderflag','no');
[traces,traceh]=getTraces(trch);
delete(trch);
end
