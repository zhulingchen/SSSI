function writeHeader ( obj )
%
%function writeHeader ( obj )
% Warning:
%   No conversions are performed, writes EBCDIC if EBCDIC, ASCII if ASCII
%

try
    %make sure we're at the start of the header in the file
    fseek(obj.fid,obj.hdroffset,'bof');

    %read text header
    fwrite(obj.fid, obj.header', '*uchar',0,obj.machineformat);

catch me
    error (me.message);
end

end

