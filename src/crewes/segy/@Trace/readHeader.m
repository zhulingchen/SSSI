function obj = readTrace ( obj )
%
%function obj = readHeader ( obj )
%
% Read the trace headers from a SEG-Y or SU file
% Returns:
%   header = a structure array containing the trace header information.
%     Each structure heading is a variable defined in headerdefininitions
%

try
    %position pointer at start of header
    fseek(obj.fid,obj.traceoffset,'bof');
    %get size of traceheader definitions
    sz=size(obj.definitions.values);
    %create structure with traceheader definintion names
    nsind=0;
    tps=cell(sz(1),1);
    obj.header.troffset=uint32([]);
    for k=1:sz(1)
        ts=obj.definitions.values{k,4};
        if strcmp(ts(1,1),'u')
            tps{k,1}=ts;
        elseif strcmp(ts(1,1),'c')
            tps{k,1}=ts;
        else
            tps{k,1}=['u',ts];
        end
        if strcmp(obj.definitions.values{k,1},'ns')
            nsind=k;
        end
    end
    header=uint32(ones(sz));
    %create variables to hold assorted things
    nsforall=[];
    val=1;
    %read trace header values
    for k=1:sz(1)
        ln=str2double(obj.definitions.values{k,3})-str2double(obj.definitions.values{k,2})+1;
        typ=obj.definitions.values{k,4};
        if findstr(typ,'32')
            ln=ln/4;
        elseif findstr(typ,'16')
            ln=ln/2;
        elseif findstr(typ,'8')
            ln=ln/1;
        end
        val=fread(obj.fid, ln ,typ);
        header(k)=val;
        
    end
    
    if nsind
        ns=header(nsind);
    end
    
    % ask user to imput number of samples if trace header does not exist
    if ~nsind || ns==0
        ns=str2double(inputdlg({['The number of samples in the trace was not found.',...
            'Please enter the number of samples in the trace.']},...
            'Number of Samples NOT FOUND',1,{'500'}));
    end
    
    obj.data=fread(obj.fid,ns,'float');
    obj.header=header;
catch me
    error (me.message);
end
end
