function expdepth
% Export the receiver coordinates, their surface elevation, the elevation of 
% the layer1-layer2 interface and the thickness of the first layer 
depth = refdata('get','depth');
recelev = refdata('get','recelev');
depthelev = recelev(2,:)-depth(2,:);
reccoord = recelev(1,:);
recnum = 1:length(reccoord);
% Ask for filename
[filename,path] = myuifile(gcf, '*.dpt', 'Export Depth', 'put');
if( filename == 0 )
   return;
end
ind = findstr(filename,'.dpt');
if(strcmp(computer,'MAC2'))
%   fullfilename=filename(1:ind(1)-1);
   fullfilename=filename;
else
%   fullfilename=[path filename(1:ind(1)-1)];
   fullfilename=[path filename];
end
fid = fopen(fullfilename, 'w');
if( fid ~= -1 )
   fprintf(fid, ' #Rec.   Coord.    Surf.     Interf.    Thick.\n');
   for i=1:length(reccoord)
      count = fprintf(fid, '%5d %9.2f %9.2f %9.2f %9.2f\n', recnum(i), ...
                      reccoord(i), recelev(2,i), depthelev(i), depth(2,i) );
   end
   fclose(fid);
else
   str = sprintf('Could not open file: %s',fullfilename);
   error(str);
end
