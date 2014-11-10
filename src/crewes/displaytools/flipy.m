function flipy
% FLIPY: script to flip the direction of the y axis
%
% just type "flipy" at the matlab prompt 
state=get(gca,'ydir');
if(strcmp(state,'normal'))
	set(gca,'ydir','reverse')
else
	set(gca,'ydir','normal')
end
