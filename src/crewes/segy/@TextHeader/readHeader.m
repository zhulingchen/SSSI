function [header format] = readHeader ( obj )
%
%function [header format] = readHeader ( obj )
%
% Read the textual file header from a SEG-Y file
% Returns:
%   thead  = 3200 element uchar 1D vector
%   format = 'ascii' or 'ebcdic', where format is set to 'ebcdic'
%            only if the header can be proven to be ebcdic
% Warning:
%   returns EBCDIC header if EBCDIC, ASCII if ASCII
%

format = 'ascii';
try
    %position pointer at start of header
    fseek(obj.fid,obj.hdroffset,'bof');

    %read text header
    header = fread(obj.fid, obj.hdrsize, '*char');
    format = obj.guessTextFormat(header);
catch me
    error (me.message);
end

end

