%% test tau-p transform code
%make a text shot record
nx=200;
nt=501;
dx=10;
dt=.002;
x=(0:nx-1)*dx;
t=((0:nt-1)*dt)';
shot=zeros(nt,nx);
v0Left=1750;
v0Right=1250;
v1=2000;
v2=3000;
xs=(nx/2)*dx;
xoffmax=max(x-xs);
t0=.1;
shot=event_dip(shot,t,x,[xoffmax/v0Left+t0 t0],[0 xs],[.1 .2]);
shot=event_dip(shot,t,x,[t0 xoffmax/v0Right+t0],[xs max(x)],[.2 .1]);
shot=event_hyp(shot,t,x,max(t)/4,xs,2*v1,.1);
shot=event_hyp(shot,t,x,max(t)/2,xs,2*v2,-.1);
[w,tw]=ricker(dt,30,.2);
shot=convz(shot,w);
xoff=x-xs;
plotimage(shot,t,xoff);
title('Shot gather to be transformed')
xlabel('Offset (m)')
v0=min([v0Left v0Right]);
pmin=-2/v0;pmax=-pmin;
dp=.25*(pmax-pmin)/nx;

[stp,tau,p]=tptran(shot,t,xoff,pmin,pmax,dp);

plotimage(stp,tau,p)
title('Tau-p transform of shot gather')
xlabel('Slowness (s/m)')

%% test inverse tau-p transform
%run the previous first
xmin=xoff(1);
xmax=xoff(end);
% dx=2*(xmax-xmin)/(length(p)-1);
[seis2,t2,x2]=itptran(stp,tau,p,xoff(1),xoff(end));

plotimage(seis2,t2,x2)
title('reconstructed shot gather')
