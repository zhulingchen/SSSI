function writetraces(fileout,traceheaders,traces,machineformat,traceoffset)
null=NaN;

try
if nargin <5
    traceoffset='3600';
end
if nargin <4
    machineformat='ieee-be';
end
szhead=size(traceheaders);
sztrac=size(traces);

if szhead(2)~=sztrac(2)
    me=MException('writetraces:SizeMismatch',...
        'There must be a traceheader entry for each trace.');
    throw(me);
end

trchnew=TraceHeader(fileout,'hdroffset',traceoffset,'machineformat',machineformat,...
    'permission','r+');

trchnew.nontypecasthdr=traceheaders;
fseek(trchnew.fid,trchnew.hdroffset,'bof');

for m=1:sztrac(2)
    
    %remove any NaNs from the trace and make sure ns is set to tracelength
    ntrace=traces(:,m);
    ntrace=ntrace(ntrace~=null);
    ns=getheadervalue(trchnew,'ns');
    ns(1,m)=length(ntrace);
    trchnew=setheadervalue(trchnew,'ns',ns);
    %write the trace header
    trchead=trchnew.nontypecasthdr(:,m);
    fwrite(trchnew.fid,trchead,'uint8');
    
    % write the trace data
    fwrite(trchnew.fid,ntrace,'float',0,trchnew.machineformat);
end

catch me
    errordlg(me.message, me.identifier);
end


end