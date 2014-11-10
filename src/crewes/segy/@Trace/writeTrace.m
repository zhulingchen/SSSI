function writeTrace( obj )
%
%function writeTrace ( obj )
% writes the traceheaders and tracedata in a trace object
%
null=NaN;
try
    %make sure we're at the start of the header in the file
    fseek(obj.traceheader.fid,obj.traceheader.hdroffset,'bof');
sz=size(obj.traceheader.nontypecasthdr);
% test to see if any nulls are in the data
if any(any(obj.tracedata.data==null))
hwait=waitbar(0,'Please Wait as Traces are Prepared for Writing');
for m=1:sz(2)
    
    %remove any null from the trace and make sure ns is set to tracelength
    ntrace=obj.tracedata.data(:,m);
    if isnan(null)
        ntrace=ntrace(~isnan(ntrace));
    else
    ntrace=ntrace(~(ntrace==null));
    end
    ns=obj.traceheader.getheadervalue('ns');
    ns(1,m)=length(ntrace);
    obj.traceheader=obj.traceheader.setheadervalue('ns',ns);
     if any(m==1:50:sz(2))
        waitbar(m/sz(2),hwait,[num2str(m),' of ',...
            num2str(sz(2)),' Traces have been Prepared for Writing. ']);
    end
end
delete(hwait)
else
    ns=size(obj.tracedata.data);
    ns=ns(1)*ones(ns(2),1);
    obj.traceheader=obj.traceheader.setheadervalue('ns',ns);
end
% swap bytes in header if nessasary
if isempty(strfind(lower(obj.traceheader.filefmt),lower(obj.traceheader.machineformat)));
    words=obj.traceheader.definitions.values(:,strcmpi(obj.traceheader.definitions.keys(),'Name'));
    for k=1:length(words)
        vari=obj.traceheader.getheadervalue(words{k},0);
        vari=swapbytes(vari);
        obj.traceheader=obj.traceheader.setheadervalue(words{k},vari,0);
    end
end

hwait=waitbar(0,'Please Wait as Traces are Written to the File');
for m=1:sz(2)
    
    %remove any null from the trace and make sure ns is set to tracelength
    ntrace=obj.tracedata.data(:,m);
    if isnan(null)
        ntrace=ntrace(~isnan(ntrace));
    else
    ntrace=ntrace(~(ntrace==null));
    end
    
    %write the trace header
    trchead=obj.traceheader.nontypecasthdr(:,m);
    fwrite(obj.traceheader.fid,trchead,'uint8',0,obj.traceheader.machineformat);
    
    % write the trace data
    fwrite(obj.traceheader.fid,ntrace,'float',0,obj.traceheader.machineformat);
    
    if any(m==1:50:sz(2))
        waitbar(m/sz(2),hwait,[num2str(m),' of ',...
            num2str(sz(2)),' Traces have been Written to the File. ']);
    end
end

catch me
    delete(hwait);
    error (me.message);
end
delete(hwait);

end

