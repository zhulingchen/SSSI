%kirchhoff time migration of thrust
clear
load thrust
%make an rms velocity model
[vrms,tv,vint]=vz2vrms(vel,z,t(2)-t(1),max(t));

%plot
figure;imagesc(x,t,vint);colorbar
title('Interval velocity in time')
xlabel('meters');ylabel('seconds')
figure;imagesc(x,t,vrms);colorbar
title('RMS Velocity in time');colorbar
xlabel('meters');ylabel('seconds')

params=nan*ones(1,12);

seismig=kirk_mig(seisf,vrms,t,x,params);

plotimage(seismig,t,x)