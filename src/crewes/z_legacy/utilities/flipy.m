state=get(gca,'ydir');
if(strcmp(state,'normal'))
	set(gca,'ydir','reverse')
else
	set(gca,'ydir','normal')
end
