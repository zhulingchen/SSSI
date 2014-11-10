function [depth] = calcdepth(plust,v1rec,v2rec)
% depth = calcdepth(plust,v1rec,v2rec)
%
% Depth calculation using the Plus Time and the layer velocities at each
% receiver location.
%
% Store all the information concerning the Depth values in a matrix of 4 rows
%
% row 1	: receiver coordinates
% row 2	: average Depth values
% row 3	: fold of the Depth values
% row 4	: standard deviation of the Depth values
critang=zeros(1,length(v1rec));
depth(1,:) = plust(1,:);
depth(3,:) = plust(3,:);
for n=1:length(v1rec)
	critang(1,n) = asin(v1rec(1,n)/v2rec(1,n));
	depth(2,n) = (plust(2,n)*v1rec(1,n))/(2*cos(critang(1,n)));
	depth(4,n) = (plust(4,n)*v1rec(1,n))/(2*cos(critang(1,n)));
end
PMTsetmenus;
