% make a synthetic reflectivity and a minimum phase wavelet
dt=.002;
tmax=2;
fdom=30;
[r,t]=reflec(tmax,dt);%synthetic reflectivity
[w,tw]=wavemin(dt,fdom,.2);%min phase wavelet
%now make stationary and nonstationary traces
s=convm(r,w);%stationary trace
Q=50;
qmat=qmatrix(Q,t,w,tw);%Q matrix for Q=50.
sQ=qmat*r;%nonstationary trace
figure
plot(t,s,t,sQ,'r')%compare the traces in the time domain
legend('Stationary',['Nonstationary, Q=' int2str(Q)])
prepfig
xlabel('Time (s)')
%now make an inverse Q matrix assuming an impulse wavlet
iqmat=invq(Q,t);%use the default tolerance
%apply to the nonstatinary seismogram
sQi=iqmat*sQ;
figure
hh=plot(t,s,t,sQ,'k',t,sQi,'r.');
set(hh(2),'color',[.5 .5 .5],'linewidth',1)
legend('Stationary seismogram','Nonstationary seismogram',...
    'Inverse Q matrix applied to nonstationary seismogram')
xlabel('Time (s)')
prepfig
%%
%compute inverse Q matrix applied to Qmat
w0=iqmat*qmat;
plotimage(w0,t,t)