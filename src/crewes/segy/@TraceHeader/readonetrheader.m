function obj = readonetrheader ( obj)
%
%function obj = readHeader ( obj )
%
% Read selected trace headers from a SEG-Y or SU file
% Input:
%   obj=TraceHeader object
%   limits= a vector containing either [min:max] where min is the minimum
%               trace number and max is the maximum trace number.
%           Or [n1 n2 n3 ...], where n1, n2 ... are trace numbers in the
%               order you which they are displayed
%
% Returns:
%   nontypecasthdr = a unsigned int8 array contining the header information
%           which can be converted usingconvertHeader or getHeader
%

try
    fseek(obj.fid,obj.hdroffset,'bof');
    val=fread(obj.fid,[obj.hdrsize,1] ,'*uint8');
    obj.nontypecasthdr=val;
    ns=obj.getheadervalue('ns');
    tracebytes=ns*obj.tracetype{1}+obj.hdrsize+obj.hdroffset;
    obj.traceoffsets=[obj.hdroffset,tracebytes];
    
    trace=fread(obj.fid,[ns*obj.tracetype{1},1] ,'*uint8');
    obj2=copy(obj);
    while ~any(trace)
        tb=obj2.traceoffsets(1);
        val=fread(obj.fid,[obj.hdrsize,1] ,'*uint8');
        obj2.nontypecasthdr=val;
        ns=obj2.getheadervalue('ns');
        obj2.traceoffsets=[tb,ns*obj.tracetype{1}+obj.hdrsize+tb];
        trace=fread(obj.fid,[ns*obj2.tracetype{1},1] ,'*uint8');
    end
    tracebuff=Trace(obj2);
    obj.tracetype{2}=tracebuff.uigettracetype(tracebuff,1);
    
    
catch me
    error (me.message);
end
end
