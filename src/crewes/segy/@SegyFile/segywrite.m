function segywrite(segyobj)
    
%set all of the fids to segyobj.fid and close the other fids;
fclose(segyobj.binaryheader.fid);segyobj.binaryheader.fid=segyobj.fid;
segyobj.binaryheader.machineformat=segyobj.machineformat;
fclose(segyobj.textheader.fid);segyobj.textheader.fid=segyobj.fid;
segyobj.textheader.machineformat=segyobj.machineformat;
for k=1:length(segyobj.extendedheader)
    fclose(segyobj.extendedheader{1,k}.fid);
    segyobj.extendedheader{1,k}.fid=segyobj.fid;
    segyobj.extendedheader{1,k}.machineformat=segyobj.machineformat;
end
fclose(segyobj.trace.traceheader.fid);segyobj.trace.traceheader.fid=segyobj.fid;
segyobj.trace.traceheader.machineformat=segyobj.machineformat;
fclose(segyobj.trace.tracedata.fid);segyobj.trace.tracedata.fid=segyobj.fid;
segyobj.trace.tracedata.machineformat=segyobj.machineformat;
    
% change last 2 lines in textheader to segy rev1 standard
c39='C39 SEG Y REV1';
c40='C40 END TEXTUAL HEADER';
% size(segyobj.textheader.header)
if length(segyobj.textheader.header) >  80
    segyobj.textheader.header(3041:3054)=c39;
    segyobj.textheader.header(3121:3142)=c40;
else
    segyobj.textheader.header(39,1:length(c39))='C39 SEG Y REV1';
    segyobj.textheader.header(40,1:length(c40))='C40 END TEXTUAL HEADER';  
end

% change values in binary header to accomodate segy revision 1 standards
    szext=length(segyobj.extendedheader); 
    segyobj.binaryheader.setheadervalue('netfh',szext);

segyobj.binaryheader=segyobj.binaryheader.setheadervalue('format',5);
segyobj.binaryheader=segyobj.binaryheader.setheadervalue('rev',bitshift(uint16(1),8));

ns=single(getheadervalue(segyobj.trace.traceheader,'ns'));
if ~isempty(ns)
    segyobj.binaryheader=segyobj.binaryheader.setheadervalue('hns',mode(ns));
    if size(unique(ns))==[1,1]
        segyobj.binaryheader=segyobj.binaryheader.setheadervalue('flen',1);
    else
        segyobj.binaryheader=segyobj.binaryheader.setheadervalue('flen',0);
    end
end

dt=single(getheadervalue(segyobj.trace.traceheader,'dt'));
if ~isempty(dt)
    segyobj.binaryheader=segyobj.binaryheader.setheadervalue('hdt',mode(dt));
end

% write the files to disk 
segyobj.textheader.writeheader();
segyobj.binaryheader.writeheader();
for k=1:length(segyobj.extendedheader)
    segyobj.extendedheader{1,k}.writeheader();
end
segyobj.trace.writetrace();

% set all of the shared fids to -1
segyobj.binaryheader.fid=-1;
segyobj.textheader.fid=-1;
segyobj.trace.traceheader.fid=-1;
segyobj.trace.tracedata.fid=-1;
for k=1:length(segyobj.extendedheader)
    segyobj.extendedheader{1,k}.fid=-1;
end




end



