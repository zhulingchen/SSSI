%% script to test cmpstack
clear
%this cell makes some finite difference shot records
nshots=5;%script will choose the shot locations to span the model
%first we build a velocity model of the burried channel structure
dx=2;
xmax=1000;
zmax=1000;
vmax=4000;
vmin=1500;
zlay=[0 200 250 350 450 500 700];
vlay=[1500 1700 2000 2100 2500 3000 3500];
zchannel=500;
wchannel=50;
hchannel=30;
vchannel=2000;
[vel,x,z,zlay,vlay]=channelmodel(dx,xmax,zmax,vmax,vmin,zchannel,wchannel,hchannel,vchannel,4,1,zlay,vlay);
figure
imagesc(x,z,vel);colorbar
title('Channel model');
% finite difference
dtstep=.0002;
dt=.004;
tmax=1;
fdom=30;
[w,tw]=wavemin(dt,fdom,.2);%source waveform

[shots,t,xshots,xrecs,shotnames]=afd_shootline(dx,vel,dt,dtstep,w,tw,tmax,nshots);

stackingchart(xshots,xrecs);

plotgathers(shots,xrecs,t,'receiver coordinate (m)','time (s)',shotnames);

linename=['channelline' int2str(nshots)];
eval(['save ' linename])

fileout='test_shootline';

SEGY_WriteShotGathers(fileout,shots,t,xshots,xrecs);
%% gain, mute, nmor and stack
linename=['channelline' int2str(nshots)];
eval(['load ' linename])
% load channelline21_Drec_2
%make an rms velocity model
velrms=vzmod2vrmsmod(vel,z,dt,tmax);
velrms2=velrms(:,1)*ones(1,length(x));%alternate model without channel
figure;
imagesc(x,t,velrms);colorbar
title('RMS Velocity model')
%loop over shot records to gain+mute+nmo correct
shots_nmo=cell(1,nshots);
xcmp=cell(1,nshots);
fold=cell(1,nshots);
%specify mute
tmute0=.1;%mute time at zero offset
tmute1=.8;%mute time at reference offset
xoffref=667;%reference offset
tmute=[tmute0 tmute1];%mute times
xmute=[0 xoffref];%mute offsets
gainpow=1;%power for gainmute
dxrecs=abs(xrecs{1}(2)-xrecs{1}(1));%receiver interval
dxcmp=dxrecs/2;%cmp interval
x0cmp=x(1);%first cmp
x1cmp=x(end);%last cmp
cmp=[dxcmp x0cmp x1cmp];
flag=1;%flag for nmor_cmp
smax=30;%stretch mute percent
xgath=[490:dxcmp:510];%location of requested output gathers
shotsg=gainmute(shots,t,xrecs,xshots,xmute,tmute,gainpow);
normopt=2;
[stack,xcmp,stackfold,gathers,xoffs]=cmpstack(shotsg,t,xrecs,xshots,cmp,velrms2,t,x,normopt,smax,xgath);
plotimage(stack,t,xcmp)
title(['Stack, normopt= ' int2str(normopt)])

%% make an exploding reflector model of the same channel
%remake the model
dx=2;
xmax=1000;
zmax=1000;
vmax=4000;
vmin=1500;
zlay=[0 200 250 350 450 500 700];
vlay=[1500 1700 2000 2100 2500 3000 3500];
zchannel=500;
wchannel=50;
hchannel=30;
vchannel=2000;
[vel2,x2,z2,zlay,vlay]=channelmodel(dx,xmax,zmax,vmax,vmin,zchannel,wchannel,hchannel,vchannel,4,1,zlay,vlay);
dtstep=.0002;
dt=.004;
tmax=1;
[zos,zosa,t]=afd_explode(dx,dtstep,dt,tmax,vel2,x2,zeros(size(x2)),w,tw,2);
plotimage(zos,t,x2);
title('Exploding reflector')