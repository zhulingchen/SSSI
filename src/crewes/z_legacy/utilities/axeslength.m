function [xlen,ylen]=axeslength(hax,unit)
% [xlen,ylen]=axeslength(hax,unit)
% [xlen,ylen]=axeslength(hax)
%
% AXESLENGTH returns x-axis and y-axis lengths for a
% constant aspect ratio. Lengths are returned in user
% or default specified units.
%
% hax	= axes handle, for example gca
% unit	= 'inches', 'centimeters', 'pixels' etc.
%=============== Default = 'centimeters' ===============
%
% xlen	= length of x-axis in units
% ylen	= length of y-axis in units
%
if nargin<2
	unit='centimeters';
end
	
	au=get(hax,'units');
	a=get(hax,'DataAspectRatio');
	a=a(1);
	set(hax,'units',unit);
	pos=get(hax,'position');
	
	if(isnan(a))
		xlen=pos(3);
		ylen=pos(4);
	elseif( a< pos(3)/pos(4) )
		xlen=pos(4)*a;
		ylen=pos(4);
	else
		ylen=pos(3)/a;
		xlen=pos(3);
	end
	
	set(hax,'units',au);
