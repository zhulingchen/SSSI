%% %do a finite-difference model of thrust
modelname='thrust model';
dx=5;
vlow=2000;vhigh=3500;
xmax=5100;zmax=2500;
[velt,x,z]=thrustmodel(dx,xmax,zmax,vhigh,vlow);
dt=.004; %temporal sample rate
dtstep=.001;
tmax=2*zmax/vlow; %maximum time
[seisfiltt,seis,t]=afd_explode(dx,dtstep,dt,tmax, ...
 		velt,x,zeros(size(x)),[5 10 40 50],0,2);
    
raymig(seisfiltt,velt,t,x,z,modelname)

%% depth migration of thrust

zcheck=0:100:2000;
[zosmig,exzos]=pspi_stack(seisfiltt,t,x,velt,x,z,[5 50],zcheck);
plotimage(zosmig,x,z);
xs=cell(size(exzos));
ts=cell(size(exzos));
titles=cell(size(exzos));
for k=1:length(exzos)
    xs{k}=(x(2)-x(1))*(0:size(exzos{k},2)-1);
    ts{k}=(t(2)-t(1))*(0:size(exzos{k},1)-1);
    titles{k}=['Extrapolated to ' int2str(zcheck(k))];
end
save thrustdata

%load the extrapolations into plotgathers
plotgathers(exzos,xs,ts,'distance (m)','time (s)',titles);