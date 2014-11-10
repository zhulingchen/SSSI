function obj=gettracepointers(obj)
% tracehead=gettracepointers(tracehead)
%   this function returns a vector
%

try
    %position pointer at start of header
    fseek(obj.fid,obj.hdroffset,'bof');
    %get size of traceheader definitions
    sz=size(obj.definitions.values);
    
    %create variables to hold assorted things
    nsforall=[];
    val=1;
    % test to find where ns (number of samples is located);
    nsind=0;
    nspos=find(strcmp(obj.definitions.values(:,strcmp(obj.definitions.keys,'Name')),'ns'));
    if nspos
        nsind=nspos;
    end
    troffset=obj.hdroffset;
    fileinfo=dir(obj.filename);
    
    %read trace header values
    m=1;  % is the number of traces that have been read in
    while(~isempty(val))
        % read in traceheader as nusigned byte intergers
        val=fread(obj.fid,[obj.hdrsize,1] ,'*uint8');
        if isempty(val)
            break
        end
        obj.nontypecasthdr=val;
        if nsind
            st=str2double(obj.definitions.values(nsind,strcmp(obj.definitions.keys,'startByte')));
            ed=str2double(obj.definitions.values(nsind,strcmp(obj.definitions.keys,'endByte')));
            typ=obj.definitions.values{nsind,strcmp(obj.definitions.keys,'Type')};
            ns=typecast(obj.nontypecasthdr(st:ed,1)',typ);
            ns=checkforrightbyteorder(ns,obj.filefmt);
        end
        
        
        % ask user to imput number of samples if trace header does not exist
        if ~nsind || ns==0 && isempty(nsforall)
            ns=str2double(inputdlg({['The number of samples in the trace was not found.',...
                'Please enter the number of samples in the trace.']},...
                'Number of Samples NOT FOUND',1,{'500'}));
            nsforall=questdlg(['You have entered ',num2str(ns),' as the number of samples in this trace.  ',...
                'Would you like to use this value for all traces in the file?'],...
                'Constant Number of Samples','Yes','No','Yes');
            if strcmp(nsforall,'No'),nsforall=[];else nsforall=ns;end;
        end
        
        if ~isempty(nsforall), ns=nsforall; end
        fseek(obj.fid,ns*4,'cof');
        troffset(1,m+1)=troffset(m)+double(240+ns*4);
        
        totcount=ceil(double((fileinfo.bytes-obj.hdroffset)/(obj.hdrsize+ns*4)));
        
        % See if additional space needs to be added to preallocate for
        % speed
        matsz=size(troffset);
        if matsz(2)<totcount;
            [user, sys] = memory;
            MemMaxBytes = sys.PhysicalMemory.Available;
            mryflag=((matsz(1)*totcount*8)<(0.1*MemMaxBytes));
            if mryflag
                troffset=[troffset,NaN(matsz(1),(totcount-matsz(2)))];
            else
                warndlg('Not Enough Memory to Read Trace Pointers')
                return
            end
        end
        
        % add one to m for the times it has gone arround the while loop
        m=m+1;
    end
    
    % remove additional entries
    troffset=troffset(:,1:m);
    
    
    % set the object header property
    obj.traceoffsets=troffset;
catch me
    error (me.message);
end
end
