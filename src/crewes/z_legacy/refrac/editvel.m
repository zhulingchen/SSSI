function editvel(edit,v,win)
% Edit the first and second layer velocities by dragging  
% (values at each receivers)  
v1rec=refdata('get','v1rec');
v2rec=refdata('get','v2rec');
recelev=refdata('get','recelev');
coord=recelev(1,:);
set(gcf,'units','pixels','position',[0 0 864 576]);
% Optional median filter applicable to velocity values
if (edit==0)
  % First layer velocity only
  if (v==1)
	v1rec=medfilt1(v1rec,win);
  end
  % Second layer velocity only
  if (v==2)
	v2rec=medfilt1(v2rec,win);
  end
  % Both first and second layer velocities
  if (v==3)
	v1rec=medfilt1(v1rec,win);
	v2rec=medfilt1(v2rec,win);
  end
end
hold on
v1handle = [];
if(length(v1rec)>1)
	v1handle = plot(coord,1000*v1rec,'g-.')
end
v2handle = [];
if(length(v2rec)>1)
	v2handle = plot(coord,1000*v2rec,'r-.')
end
xlabel('Coordinate (m)')
ylabel('Velocity (m/s)')
title('Velocity model (first layer=green; second layer=red)')
% Call the edit function
editveldepth('init',[v1handle v2handle],'vel');
