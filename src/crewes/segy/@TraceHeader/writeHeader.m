function writeHeader ( obj )
%
%function writeHeader ( obj )
% Warning:
%   No byte-swapping is performed
%

try
    %make sure we're at the start of the header in the file
    fseek(obj.fid,obj.hdroffset,'bof');

    %read text header
    fwrite(obj.fid, obj.hdrsize, '*uchar');
catch me
    error (me.message);
end

end

