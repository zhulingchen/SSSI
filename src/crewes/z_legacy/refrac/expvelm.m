function expvelm
recelev = refdata('get','recelev');
v1rec = refdata('get','v1rec');
v2rec = refdata('get','v2rec');
reccoord = recelev(1,:);
recnums = 1:length(reccoord);
% Ask for filename
[filename,path] = myuifile(gcf, '*.mat', 'Export Velocity', 'put');
if( filename == 0 )
   return;
end
ind = findstr(filename,'.mat');
if(strcmp(computer,'MAC2'))
%   fullfilename=filename(1:ind(1)-1);
   fullfilename=filename;
else
%   fullfilename=[path filename(1:ind(1)-1)];
   fullfilename=[path filename];
end
fid = fopen(fullfilename, 'w');
count = fprintf('%d %f %f %f', recnums, reccoord, v1rec, v2rec);
