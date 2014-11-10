%exploding reflector model of thrust

dt=.004;
dtstep=.001;
tmax=2.0;
dx=5;

[vel,x,z]=thrustmodel(dx);

[seisf,seis,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),[5 10 40 50],0,2);

figure
imagesc(x,z,vel);
h=colorbar;
set(get(h,'ylabel'),'string','m/s')
title('velocity model')
plotimage(seisf,t,x)
