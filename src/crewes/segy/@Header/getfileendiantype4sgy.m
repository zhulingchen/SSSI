function endi=getfileendiantype4sgy(obj,offset)


% read in indication variable if using sgy standard bytes 3225-3226 should
% be 1-8
fid=obj.fid;
fseek(fid,offset,'bof');
num=fread(fid,1,'*int16');

% get byte order of operating system
    [ab, ac, e] = computer;
% set change endian to false
 changeendi=false;
 
% set changeendi to true if the value of the indication variable is greater
% than 256

if abs(num)>=256
    changeendi=true;
end
 
% if changeendi is true make endi the opposite of the current computer else
% if changeendi is false make endi the same as the current computer
if changeendi
    if (strcmp(e,'L'))
        endi='B';
    else
      endi='L';
    end
else
    endi=e;
end

end