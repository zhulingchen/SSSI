function expstat
% Export the receiver and the shot coordinates with their corresponding
% weathering, elevation and total static corrections
recelev = refdata('get','recelev');
recstat = refdata('get','recstat');
shotcoord = refdata('get','shotcoord');
shotstat = refdata('get','shotstat');
reccoord = recelev(1,:);
recnum = 1:length(reccoord);
shotnum = 1:length(shotcoord);
% Ask for filename
[filename,path] = myuifile(gcf, '*.sta', 'Export Static', 'put');
if( filename == 0 )
   return;
end
ind = findstr(filename,'.sta');
if(strcmp(computer,'MAC2'))
%   fullfilename=filename(1:ind(1)-1);
   fullfilename=filename;
else
%   fullfilename=[path filename(1:ind(1)-1)];
   fullfilename=[path filename];
end
fid = fopen(fullfilename, 'w');
if( fid ~= -1 )
   fprintf(fid, ' #Rec.     Coord.     Weath.      Elev.      Total\n');
   for i=1:length(reccoord)
      count = fprintf(fid, '%5d %11.2f %11.5f %11.5f %11.5f\n', ...
           recnum(i), reccoord(i), recstat(1,i), recstat(2,i), recstat(3,i) );
   end
   fprintf(fid, ' #Shot     Coord.     Weath.      Elev.      Total\n');
   for i=1:length(shotcoord)
      count = fprintf(fid, '%5d %11.2f %11.5f %11.5f %11.5f\n', recnum(i), ...
                      shotcoord(i), shotstat(1,i), shotstat(2,i), shotstat(3,i) );
   end
   fclose(fid);
else
   str = sprintf('Could not open file: %s',fullfilename);
   error(str);
end
