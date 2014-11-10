function obj = guessByteOrder(obj)
%
% Guess if data in SEG-Y file are big or little endian
% Theory:
%   The data format code in the binary file header
%   will be out of range [1,8] if read using
%   the wrong byte order.
% Returns:
%   'ieee-be' (SEG-Y standard) unless the data can be proven to be
%   little-endian (PC byte order), in which case 'ieee-le' is returned
%

%Assume standard segy
obj.machineformat = 'ieee-be';

%Test for PC byte order
try
    fseek(obj.fid,3224,'bof');
    if(fread(obj.fid, 1, 'int16', 0, 'ieee-le') <255)
        obj.machineformat = 'ieee-le';
    end
catch me
    error(me.message);
end

end

