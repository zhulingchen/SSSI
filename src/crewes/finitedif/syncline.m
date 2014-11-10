%%
global SCALE_OPT
SCALE_OPT=2;
% buried synlcine (focal point below the surface)
%see help for synclinemodel to understand the next few parameters
dx=5;xmax=2000;
zmax=1000;
vhigh=4000;
vlow=2000;
zsyncline=300;
zfocal=100;
radius=500; %radius+zfocal should exceed zsyncline
[vel,x,z]=synclinemodel(dx,xmax,zmax,vhigh,vlow,zsyncline,zfocal,radius);
figure
imagesc(x,z,vel);colorbar
title(['Syncline model zfocal= ' num2str(zfocal) 'm'])
% prepare for finite difference modelling. See help for afd_explode
dt=.004;
dtstep=.0005;
tmax=1.5;
%wavelet
[w,tw]=ricker(dt,30,.2);
%
[seis1,seis_unfiltered,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
% [seis,t]=afd_explode_alt(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
plotimage(seis1,t,x)
title(['Syncline with focal depth at z= ' num2str(zfocal) 'm'])
%%
% focal point above surface
global SCALE_OPT
SCALE_OPT=2;
dx=5;xmax=2000;
zmax=1000;
vhigh=4000;
vlow=2000;
zsyncline=300;
zfocal=-100;
radius=500; %radius+zfocal should exceed zsyncline
[vel,x,z]=synclinemodel(dx,xmax,zmax,vhigh,vlow,zsyncline,zfocal,radius);
figure
imagesc(x,z,vel);colorbar
title(['Syncline model zfocal= ' num2str(zfocal) 'm'])
% prepare for finite difference modelling
dt=.004;
dtstep=.0005;
tmax=1.5;
%wavelet
[w,tw]=ricker(dt,30,.2);
%
[seis2,seis_unfiltered,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
% [seis,t]=afd_explode_alt(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
plotimage(seis2,t,x)
title(['Syncline with focal depth at z= ' num2str(zfocal) 'm'])
%%
% buried sinusoidal syncline
global SCALE_OPT
SCALE_OPT=2;
dx=5;xmax=2000;
x=0:dx:xmax;
zmax=1000;
z=0:dx:zmax;
vhigh=4000;
vlow=2000;
wavelength=.5*xmax;
zsyncline=500;%mean depth of the sinusoid
zamp=200;%amplitude of the sinusoid.
ztop=zamp*cos((x-mean(x))*2*pi/(wavelength))+zsyncline;
xpoly=[x fliplr(x)];
zpoly=[zeros(size(x)) ztop];
vel=vhigh*ones(length(z),length(x));
vel=afd_vmodel(dx,vel,vlow,xpoly,zpoly);
figure
imagesc(x,z,vel);colorbar
title(['Co-Syncline model'])
% prepare for finite difference modelling
dt=.004;
dtstep=.0005;
tmax=1.5;
%wavelet
[w,tw]=ricker(dt,30,.2);
%
[seis3,seis_unfiltered,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
% [seis,t]=afd_explode_alt(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
plotimage(seis3,t,x)
title(['Co-Syncline model'])
%%
% un-buried sinusoidal syncline
global SCALE_OPT
SCALE_OPT=2;
dx=5;xmax=2000;
x=0:dx:xmax;
zmax=1000;
z=0:dx:zmax;
vhigh=4000;
vlow=2000;
wavelength=.5*xmax;
zsyncline=150;%mean depth of the sinusoid
zamp=100;%amplitude of the sinusoid.
ztop=zamp*cos((x-mean(x))*2*pi/(wavelength))+zsyncline;
xpoly=[x fliplr(x)];
zpoly=[zeros(size(x)) ztop];
vel=vhigh*ones(length(z),length(x));
vel=afd_vmodel(dx,vel,vlow,xpoly,zpoly);
figure
imagesc(x,z,vel);colorbar
title(['Co-Syncline model'])
% prepare for finite difference modelling
dt=.004;
dtstep=.0005;
tmax=1.5;
%wavelet
[w,tw]=ricker(dt,30,.2);
%
[seis3,seis_unfiltered,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
% [seis,t]=afd_explode_alt(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),w,tw,2);
plotimage(seis3,t,x)
title(['Co-Syncline model'])