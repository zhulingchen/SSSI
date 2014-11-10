function obj = readHeader ( obj )
%
%function obj = readHeader (obj)
%
%   Reads header from SEG-D file as type 'uint8'
%

obj.header    = fread( obj.fid, obj.headerSize, 'uint8=>uint8');
obj.hexHeader = sprintf('%02X',obj.header);

end



