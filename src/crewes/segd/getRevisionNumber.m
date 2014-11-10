function [rmaj rmin] = getRevisionNumber(fid)
%
%function [rmaj rmin] = getRevisionNumber (fid)
%
% The SEG-D revision number, for revisions 1+, should be at bytes 11 and 12
% of General header block 2, or bytes 43 and 44 in the file
%
% fid = fopen( ... ) needs to be used before calling this function
%

%get current position in file
filepos = ftell(fid);

%move to where revision number is located
fseek(fid,42,'bof');

%get revision number
rmaj = fread(fid,1,'uint8=>uint8');
rmin = fread(fid,1,'uint8=>uint8');

%reset file position
fseek(fid,filepos,'bof');

end



