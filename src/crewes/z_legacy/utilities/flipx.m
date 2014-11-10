state=get(gca,'xdir');
if(strcmp(state,'normal'))
	set(gca,'xdir','reverse')
else
	set(gca,'xdir','normal')
end
