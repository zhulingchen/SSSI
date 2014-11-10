function tracetype=gettracetype(obj,offset)

tracetype={4,'unknown'};

if nargin<2
    offset=3224;
end
% get tracetype flag from binary header
%for segy standard offset=3224;
fseek(obj.fid,offset,'bof');
num=fread(obj.fid,1,'*int16');
num=checkforrightbyteorder(num,obj.filefmt);

if     num==1
    tracetype={4,'unknown'};
elseif num==2
    tracetype={4,'int32'};
elseif num==3
    tracetype={2,'int16'};
elseif num==4
    tracetype={4,'float32'};
elseif num==5
    tracetype={4,'float32'};
elseif num==8
    tracetype={1,'int8'};
else
    tracetype={4,'unknown'};
end



end