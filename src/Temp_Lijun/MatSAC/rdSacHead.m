function hd = rdSacHead(sacFile)
% read header of SAC format data

aa=fopen(sacFile,'r');
if aa<0
   sacFile
   return;
end
part1 = fread(aa,70,'float');
part2 = fread(aa,40,'int');
fclose(aa);

hd.delta = part1(1);
hd.b	 = part1(6);
hd.stla	 = part1(32);
hd.stlo	 = part1(33);
hd.evla	 = part1(36);
hd.evlo	 = part1(37)
hd.dist	 = part1(51);
hd.az	 = part1(52);
hd.baz	 = part1(53);
hd.gcarc = part1(54);
hd.npts	 = part2(10);
