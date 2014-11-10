% Function: readPMdb
%
% Usage:  [num value] = readPMdb('filename.a_db');
%
% Reads a ProMAX ASCII database output file, returning two vectors
% containing the values of that database order.  For example, if
% you export the receiver X coordinate database, the first vector
% (num) will be the receiver number, and the 2nd vector (value) is
% the receiver coordinate.
function [a, b] = readPMdb(filename)
if( nargin < 1 )
   error('readPMdb: No file name specified.');
end
fid = fopen(filename,'r');
if( fid == -1 )
   errs = sprintf('readPMdb: Could not open file %s\n',filename);
   error(errs);
end
% A ProMAX database export file looks like:
%
%
%
%  ASCII database file write for Area="blackfoot", Line="10hzp"
%
%
%  Value list for Database Order SRF
% Surface location X coordinates                                                  
%>SRFGEOMETRYX_COORD 
%        SRF   X_COORD 
%<       101| 348991.22 |
%<       102| 348977.50 |
% Read and skip over the the first 10 lines
for i=1:10
   tmp = fgetl(fid);
end
% Next, read each line, until we don't successfully read 2 numbers
count = 2;
i = 1;
[values count] = fscanf(fid, '< %d| %f |');
pairs = reshape(values,2,count/2)';
a = pairs(:,1);
b = pairs(:,2);
fprintf(1,'readPMdb: read %d data pairs from %s\n',count/2,filename);
fclose(fid);
