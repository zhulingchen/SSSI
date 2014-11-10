function obj = readlimitedHeader ( obj, limits)
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
    troffset=obj.traceoffsets;
    val=1;
    sz=size(limits);
    % check for enough memory
    [user, sys] = memory;
            MemMaxBytes = sys.PhysicalMemory.Available;
            mryflag=((obj.hdrsize*sz(2))<(0.1*MemMaxBytes));
    if mryflag
        nontypecastheader=uint8(zeros(obj.hdrsize, sz(2)));
    else
        warndlg('Not Enough Memory to Read Trace Headers')
        return
    end
    hwait=waitbar(0/length(limits),'Please Wait as Trace Headers Load');
    
    
    %read trace header values
    for m=1:sz(2)
        % read in traceheader as nusigned byte intergers
        fseek(obj.fid,troffset(limits(m)),'bof');
        val=fread(obj.fid,[obj.hdrsize,1] ,'*uint8');
        if isempty(val)
            break
        end
        nontypecastheader(:,m)=val;
        
    end
    
    % set the object header property
    obj.nontypecasthdr=nontypecastheader;
catch me
    error (me.message);
end
delete(hwait);
end
