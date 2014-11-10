function flipx
% FLIPX: script to flip the direction of the x axis
%
% just type "flipx" at the matlab prompt 
state=get(gca,'xdir');
if(strcmp(state,'normal'))
	set(gca,'xdir','reverse')
else
	set(gca,'xdir','normal')
end
