function obj = convertHeader ( obj )
%
%function obj = convertHeader ( obj )
%
% Refomats the header information as defined by HeaderDefinitions
% Returns:
%   header = a structure array containing the trace header information.
%     Each structure heading is a variable defined in headerdefininitions
%

try
    %size of definitions
    sz=size(obj.definitions.values);
    numoftr=size(obj.nontypecasthdr);
    hwait=waitbar(0,'Please Wait as Trace Headers are Converted');
    % loop to typecast unformatted traceheaders
    
    for k=1:sz(1)
        st=str2double(obj.definitions.values(k,strcmp(obj.definitions.keys,'startByte')));
        ed=str2double(obj.definitions.values(k,strcmp(obj.definitions.keys,'endByte')));
        typ=obj.definitions.values{k,strcmp(obj.definitions.keys,'Type')};
        ibmflag=false;
        if strfind(typ,'ieee')
            typ='float32';
        elseif strfind(typ,'ibm')
            typ='uint32';
            ibmflag=true;
        end
        val=eval([typ,'(zeros(1,','numoftr(2)))']);    
        for m=1:numoftr(2)
            val(1,m)=typecast(obj.nontypecasthdr(st:ed,m)',typ);
        end
        val=checkforrightbyteorder(val,obj.filefmt);
        if ibmflag
            val=ibm2ieee(val);
        end
        header.(obj.definitions.values{k,1})=val;
        %Adjust waitbar
        waitbar(k/sz(1));
        
    end

% set the object header property
obj.header=header;

catch me
    error (me.message);
end
delete(hwait);
end
