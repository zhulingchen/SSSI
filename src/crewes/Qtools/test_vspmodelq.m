%% Cell 1, a very simple model
%make a very simple model
z=(0:5:2000)';
vp=2000*ones(size(z));
%install an interfaces at 250m and at 750m
i1=100/5+1;
i2=750/5+1;
vp(i1:i2)=2000;
vp(i2:end)=3000;
rho=2000*ones(size(z));
rho(i1:i2)=2000;
rho(i2:end)=2200;
Q=50*ones(size(z));
Q(i1:i2)=50;
Q(i2:end)=100;

figure
plot(vp,z,rho,z,100*Q,z);flipy
title('Simple model')
ylabel('Depth (m)');
xlabel('Velocity or Density (MKS units)')
legend('Velocity','Density','100*Q');
axis([1500 10000 0 2000])
prepfig

dt=.002;
fdom=30;
% [w,tw]=ricker(dt,30,.2);
[w,tw]=wavemin(dt,fdom,.2);
figure
plot(tw,w)
title('Wavelet')
xlabel('time (s)')
prepfig
zr=0:10:1000;% receiver depths

%compute a time depth curve from the well velocities
tz=vint2t(vp,z);
%interpolate times at each receiver depth
trec=interp1(z,tz,zr);
%compute drift times
td=tdrift(Q,z,fdom,vp);
%interpolate at the receiver depth
tdr=interp1(z,td,zr);

%run the model
tmax=1.5;
rflag=0;
fpress=0;
fmult=0;
[vspq,tq,upq,downq]=vspmodelq(vp,rho,Q,z,w,tw,tmax,zr,12500,rflag,fpress,fmult);

plotimage(vspq,tq,zr);
title(['Total field, rflag=' num2str(rflag)] );
xlabel('receiver depth');
ylabel('time');
plotimage(upq,tq);
title(['Upgoing field, rflag=' num2str(rflag)] );
xlabel('receiver depth');
ylabel('time');
plotimage(downq,tq,zr);
h1=line(zr,trec,'color','r');
h2=line(zr,trec+tdr,'color','g');
title(['Downgoing field, rflag=' num2str(rflag)]);
xlabel('receiver depth');
ylabel('time');
legend([h1 h2],'time at well velocity','time at seismic velocity');
figure
trace=vspq(:,end);
amp=max(abs(trace))/3;
ind=near(tq,trec(end-1));
p1=zeros(size(tq));
p1(ind)=amp;
ind=near(tq,trec(end-1)+tdr(end-1));
p2=zeros(size(tq));
p2(ind)=amp;
plot(tq,trace,tq,p1,'r',tq,p2,'k')
xlabel('Time (s)')
legend('Deepest receiver','Well velocity time','Seismic velocity time')


%% Cell 2: Read well 1409 and block it in prep for a more complex model
%load a well log and block it
wellname='1409.las';%las file to read
s=which('vspmodelq');
if(isempty(s))
    error('Well file for well 1409 not found, you need to load and install the CREWES toolbox')
end
ind = strfind(s,'vspmodelq');
filename=[s(1:ind-1) '1409.las'];
disp(['Well 1409 loaded from ' filename])
dzblk=1;%blocking size
dzout=1;%sample size
vp0=1600;
vs0=900;
rho0=1800;

[vp,vs,rho,z]=blocklogs(filename,dzblk,dzout,vp0,vs0,rho0);

%invent a Q
Qmin=20;Qmax=500;
vp1=1500;vp2=4500;
rho1=1800;rho2=3000;
[Q,Qrand]=fakeq(vp,rho,Qmin,Qmax,2,vp1,vp2,rho1,rho2,1,1);

figure
plot(vp,z,rho,z,10*Q,z);flipy
title('Well model')
ylabel('Depth (m)');
xlabel('Velocity or Density (MKS units) and Q')
legend('Velocity','Density','10*Q');
prepfig
%% Cell 3: create the VSP on the blocked logs model from Cell 2
%run the well model from the previous cell

%define wavelet and receivers
dt=.002;
fdom=30;
% [w,tw]=ricker(dt,30,.2);
[w,tw]=wavemin(dt,fdom,.2);

figure
plot(tw,w)
title('Wavelet')
xlabel('time (s)')
prepfig

tmax=2.0;
zr=0:5:1400;
rflag=1;
fpress=0;
fmult=1;
[vspq,tq,upq,downq]=vspmodelq(vp,rho,Q,z,w,tw,tmax,zr,12500,rflag,fpress,fmult);

%compute a time depth curve from the well velocities
tz=vint2t(vp,z);
%interpolate times at each receiver depth
trec=interp1(z,tz,zr);
%compute drift times
td=tdrift(Q,z,fdom,vp);
%interpolate at the receiver depth
tdr=interp1(z,td,zr);

plotimage(vspq',zr,tq);
% h1=line(zr,tr,'color','r');
% h2=line(zr,tr+tdr,'color','g');
title('Total field');
ylabel('receiver depth');
xlabel('time');


% plotimage(vspq,tq,zr);
% % h1=line(zr,tr,'color','r');
% % h2=line(zr,tr+tdr,'color','g');
% title('Total field');
% xlabel('receiver depth');
% ylabel('time');
% legend([h1 h2],'time at well velocity','time at seismic velocity');
% figure
% trace=vspq(:,end);
% amp=max(abs(trace))/3;
% ind=near(tq,trec(end-1));
% p1=zeros(size(tq));
% p1(ind)=amp;
% ind=near(tq,trec(end-1)+tdr(end-1));
% p2=zeros(size(tq));
% p2(ind)=amp;
% plot(tq,trace,tq,p1,'r',tq,p2,'k')
% xlabel('Time (s)')
% legend('Deepest receiver','Well velocity time','Seismic velocity time')
plotimage(upq',zr,tq);
title('Upgoing field');
ylabel('receiver depth');
xlabel('time');
plotimage(downq',zr,tq);
h1=line(trec,zr,'color','r');
h2=line(trec+tdr,zr,'color','g');
title('Downgoing field');
ylabel('receiver depth');
xlabel('time');
legend([h1 h2],'time at well velocity','time at seismic velocity');
figure
trace=downq(:,end);
amp=max(abs(trace))/3;
ind=near(tq,trec(end-1));
p1=zeros(size(tq));
p1(ind)=amp;
ind=near(tq,trec(end-1)+tdr(end-1));
p2=zeros(size(tq));
p2(ind)=amp;
plot(tq,trace,tq,p1,'r',tq,p2,'k')
title('Downgoing field at deepest receiver') 
xlabel('Time (s)')
legend('Deepest receiver','Well velocity time','Seismic velocity time')

%save the results
save synthetic_vsp_1409