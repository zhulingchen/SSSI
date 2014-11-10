function traceheaderposition=testforextendedheaders(obj)
fid=obj.fid;
filein=obj.filename;
trstat=obj.hdroffset;

if nargin<2
    trstat=3600;
end



fseek(fid,3504,'bof');
numofext=fread(fid,1,'uint16');
if numofext
    if numofext==-1
        numofext=500;
    end
    for n=1:numofext
        extheadtmp=TextHeader(filein,'txthoffset',num2str(trstat));
        response=questdlg(extheadtmp.header,'Is This An Extended Header?','yes','No','No');
        if strcmp(response,'yes')
            trstat=trstat+3200;
        else
            break;
        end
    end
    
end
traceheaderposition=trstat;
end
