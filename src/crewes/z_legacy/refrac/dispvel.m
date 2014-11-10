function dispvel
% Display of the first and second layer velocity for each receiver
f = gcf;
v1rec=refdata('get','v1rec');
v2rec=refdata('get','v2rec');
recelev=refdata('get','recelev');
coord=recelev(1,:);
figure
hold on
if(length(v1rec)>1)
	plot(coord,1000*v1rec,'g+')
end
if(length(v2rec)>1)
	plot(coord,1000*v2rec,'rx')
end
xlabel('Coordinate (m)')
ylabel('Velocity (m/s)')
title('Velocity model (first layer= +green; second layer= xred)')
set(gcf,'units','pixels','position',[0 0 864 576],'menubar','none');
figure(f);
