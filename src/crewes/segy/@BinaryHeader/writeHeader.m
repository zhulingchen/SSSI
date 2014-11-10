function writeHeader ( obj )
%
%function writeHeader ( obj )
% Warning:
%   No byte-swapping is performed
%

try
    %make sure we're at the start of the header in the file
    fseek(obj.fid,obj.hdroffset,'bof');

    % swap bytes if nessasary
if isempty(strfind(lower(obj.filefmt),lower(obj.machineformat)));
    words=obj.definitions.values(:,strcmpi(obj.definitions.keys(),'Name'));
    for k=1:length(words)
        vari=obj.getheadervalue(words{k});
        vari=swapbytes(vari);
        obj=obj.setheadervalue(words{k},vari);
    end
end
    
    
    %write the binary header
    fwrite(obj.fid,obj.nontypecasthdr, 'uint8',0,obj.machineformat);

catch me
    error (me.message);
end

end

