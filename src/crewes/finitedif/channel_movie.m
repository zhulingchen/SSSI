% CHANNEL: model a channel beneath a few layers
%
% low velocity channel beneath a v(z) medium
% Just run the script

width=100;

dx=5;xmax=2500;zmax=1000;%grid size, max line length, max depth
x=0:dx:xmax;z=0:dx:zmax; % x and z coordinate vector
vhigh=3000;vlow=2000;vrange=vhigh-vlow; % high and low velocities
vel=vlow*ones(length(z),length(x));%initialize velocity matrix
z1=100;z2=200;v1=vlow+vrange/5;%first layer
xpoly=[-dx xmax+dx xmax+dx -dx];zpoly=[z1 z1 z2 z2];
vel=afd_vmodel(dx,vel,v1,xpoly,zpoly);%install first layer
z3=271;v2=vlow+2*vrange/5;zpoly=[z2 z2 z3 z3];%second layer
vel=afd_vmodel(dx,vel,v2,xpoly,zpoly);%install second layer
z4=398;v3=vlow+pi*vrange/5;zpoly=[z3 z3 z4 z4];%third layer
vel=afd_vmodel(dx,vel,v3,xpoly,zpoly);%install third layer
zpoly=[z4 z4 zmax+dx zmax+dx];%last layer
vel=afd_vmodel(dx,vel,vhigh,xpoly,zpoly);%install last layer
thk=50;vch=vlow+vrange/6;%channel
xpoly=[xmax/2-width/2 xmax/2+width/2 xmax/2+width/2 xmax/2-width/2];
zpoly=[z4 z4 z4+thk z4+thk];
vel=afd_vmodel(dx,vel,vch,xpoly,zpoly);%install channel
plotimage(vel-.5*(vhigh+vlow),z,x);%plot the velocity model

%do a finite-difference exploading reflector model
dt=.004; %temporal sample rate
dtstep=.001; %modelling step size
tmax=2*zmax/vlow; %maximum time
snap1=zeros(size(vel));
snap2=snap1;
snap2(1,75)=1;
M=afd_movie(dx,dtstep,tmax,vel,snap1,snap2,2,40);

%a shot record left of the channel
%snap1=zeros(size(vel));
%snap2=snap1;
%snap2(1,75)=1;
%[shotf,shot,tshot]=afd_shotrec(dx,dtstep,-dt,tmax,...
%	vel,snap1,snap2,x,zeros(size(x)),[10 15 40 50],0,1);

%plot the seismogram
%plotimage(shotf,tshot,x)


%a shot record over the channel
%snap1=zeros(size(vel));
%snap2=snap1;
%snap2(1,125)=1;
%[shotf2,shot2,tshot]=afd_shotrec(dx,dtstep,-dt,tmax,...
%	vel,snap1,snap2,x,zeros(size(x)),[10 15 40 50],0,1);


%plot the seismogram
%plotimage(shotf2,tshot,x)
