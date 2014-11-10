function [data, hd] = rdSac(sacFile)
% read SAC format data

aa=fopen(sacFile,'r');
if aa<0
  sacFile
  return;
end
hd = fread(aa,70,'single');
hd(71:158) = fread(aa,88,'int');
data = fread(aa,inf,'single');
fclose(aa);
