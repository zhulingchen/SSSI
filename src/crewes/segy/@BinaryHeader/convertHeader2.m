function obj=convertHeader2(obj)

%
%function obj = convertHeader ( obj )
%
% Formats the raw uint8 header according to the definitions set by
% HeaderDefinitions
%

try
    % convert to formatted binary header according to definitions
    sz=size(obj.definitions.values);
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
        val=typecast(obj.nontypecasthdr(st:ed),typ);
        val=checkforrightbyteorder(val,obj.filefmt);
        if ibmflag
            val=ibm2ieee(val);
        end
        obj.header.(obj.definitions.values{k,1})=val;
    end
    
catch me
    error (me.message);
end

end