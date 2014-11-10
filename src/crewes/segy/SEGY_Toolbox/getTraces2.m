function [traces, traceh]=getTraces(tracehead,varargin)
traces=[];
traceh=[];

try
    % test if tracehead is a TraceHeader Object
    if ~isa(tracehead,'TraceHeader');
        me=MException('getTraces:InvalidInputType',...
            ['tracehead MUST be a TraceHeader Object',...
            'Please Run getTraceHeader First']);
        throw(me)
     end
    % add the end of file placement to troffset
    troffset=tracehead.traceoffsets;
    fseek(tracehead.fid,0,'eof');
    endoffile=ftell(tracehead.fid);
    if troffset(end)~=endoffile
        troffset=[troffset,endoffile];
    end
    % find out what type of floating point 
    if strcmp(tracehead.tracetype{2},'unknown')
            tracehead.tracetype{2}=viewfloatingpointtypes(tracehead);
    end
    % test what fields the user has inputed and get the values
    
    if nargin>1
        if isodd(length(varargin))
            me=MException('getTraces:insufficentinput',...
                'Input must be entered in pairs.  ');
            throw(me);
        end
        name=cell(length(varargin)/2,1);
        value=cell(length(varargin)/2,1);
        m=1;
        for i = 1:2:length(varargin)
            name{m,1} = lower(varargin{i});
            if ~ischar(name{m,1})
               me=MException('getTraces:InputIsWrongDataType',...
                'Searching parameters must be strings.');
            throw(me); 
            end
            value{m,1}=varargin{i+1};
            if ~isnumeric(value{m,1})
               me=MException('getTraces:InputIsWrongDataType',...
                'Searching parameters must be accompanied by numeric values.');
            throw(me); 
            end
            
            m=m+1;
        end
    end
    
    fldnames=fieldnames(tracehead.header);
    masterind=[];
    %find traces that satisfy name, values
    if nargin>1
        numofth=length(getfield(tracehead.header,fldnames{2,1}));
        masterind= true(length(name),numofth);
        for i=1:length(name)
            if strcmp(name{i,1},'tracenumber')
                masterind(i,:)= false(1,numofth);
                masterind(i,value{i,1})=1;
            else
                defpos=find(strcmp(fldnames,name{i,1}));
                if defpos
                    vals=value{i,1};
                    trval=getfield(tracehead.header,fldnames{defpos,1});
                    ind= false(length(vals),numofth);
                    for n=1:length(vals)
                        ind(n,:)=trval==vals(n);
                    end
                    if n>1
                        ind=logical(sum(ind));ind(ind)=1;
                    end
                    masterind(i,:)=ind;
                end
            end
        end
        tracenum=ones(1,length(masterind));
        szmast=size(masterind);
        for n=1:szmast(1)
            tracenum=tracenum.*masterind(n,:);
        end
        tracenum=find(tracenum==1);
        
    else
        tracenum=(1:length(troffset)-1);
    end
    
    % Preallocate traces and traceh matrix
    tracelength=(troffset(2)-troffset(1)-tracehead.hdrsize)/tracehead.tracetype{1};
    traces=NaN(tracelength,length(tracenum));
    traceh=uint32(ones(length(tracehead.definitions.values),length(tracenum)));
    hwait=waitbar(0,'Please Wait as Traces are Loaded');
    
    
    for j=1:length(tracenum)
        % Read in the traces
        k=tracenum(j);
        fseek(tracehead.fid,troffset(k)+tracehead.hdrsize,'bof');
        tracelength=(troffset(k+1)-troffset(k)-tracehead.hdrsize)/tracehead.tracetype{1};
        if ~fpequal(tracelength,length(traces(:,j)),100)
            appnan=NaN((length(tracelength)-length(traces(:,j))),length(traces(1,:)));
            traces=[traces;appnan];
        end
        if strcmp(tracehead.tracetype{2},'ibm')
            trace1=fread(tracehead.fid,tracelength,'*uint32');
            trace1=checkforrightbyteorder(trace1,tracehead.filefmt);
             trace1 = ibm2ieee (trace1);
             traces(:,j)=trace1;
        else
        trace1=fread(tracehead.fid,tracelength,tracehead.tracetype{2},0,tracehead.machineformat);
        %trace1=checkforrightbyteorder(trace1,tracehead.filefmt);
        traces(:,j)=trace1;
        end
        % Extract Trace Header Data from TraceHeader Object
        
        for n=1:length(tracehead.definitions.values)
            headdat=getfield(tracehead.header,tracehead.definitions.values{n,1});
            traceh(n,j)=headdat(1,k);
        end
        
        waitbar(j/length(tracenum))
    end
    
   
    
catch me
    error(me.message)
end

delete(hwait);
end