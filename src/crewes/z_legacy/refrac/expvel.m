function expvel
% export the receiver coordinates, the first and the second layer velocities 
% at the receiver location 
recelev = refdata('get','recelev');
v1rec = 1000 * refdata('get','v1rec');
v2rec = 1000 * refdata('get','v2rec');
reccoord = recelev(1,:);
recnum = 1:length(reccoord);
% Ask for filename
[filename,path] = myuifile(gcf, '*.vel', 'Export Velocity', 'put');
if( filename == 0 )
   return;
end
ind = findstr(filename,'.vel');
if(strcmp(computer,'MAC2'))
%   fullfilename=filename(1:ind(1)-1);
   fullfilename=filename;
else
%   fullfilename=[path filename(1:ind(1)-1)];
   fullfilename=[path filename];
end
fid = fopen(fullfilename, 'w');
if( fid ~= -1 )
   fprintf(fid, '#Rec.   Coord.   V1(m/s)   V2(m/s)\n');
   for i=1:length(reccoord)
      count = fprintf(fid, '%5d %11.5f %11.5f %11.5f\n', recnum(i), ...
                      reccoord(i), v1rec(i), v2rec(i) );
   end
   fclose(fid);
else
   str = sprintf('Could not open file: %s',fullfilename);
   error(str);
end
