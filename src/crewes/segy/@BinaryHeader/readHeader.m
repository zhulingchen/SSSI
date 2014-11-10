function obj = readHeader ( obj )
%
%function obj = readHeader ( obj )
%
% Read the binary file header from a SEG-Y file
% Returns:
%   header = 400 element uchar 1D vector
%   format = 'ieee-le' or 'ieee-be'
%

try
    %position pointer at start of header
    fseek(obj.fid,obj.hdroffset,'bof');
    
    
    %read uint8 header
    obj.nontypecasthdr=fread(obj.fid,obj.hdrsize,'*uint8');
catch me
    error (me.message);
end

end
