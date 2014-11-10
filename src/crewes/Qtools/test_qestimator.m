%% test qestimator
%a single simulation using only a Q matrix (so transmission only)
dt=.002;
tmax=2;
Qnot=50;
t=0:dt:tmax;
fdom=30;
[w,tw]=wavemin(dt,fdom,.4,5);
Q=Qnot*ones(size(t));
flag=1;
qmat=qmatrix(Q,t,w,tw,flag);

t1=.5;
t2=1.0;
i1=round(t1/dt)+1;
i2=round(t2/dt)+1;
twin=.3;
f1=10;
f2=50;
wintype='gaus';
spectype='multi';
pflag=1;
method='specrat';
Qint_true=Qnot;
Qint=qestimator(qmat(:,i1),qmat(:,i2),t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
prepfig


t1=.5;
t2=1.0;
i1=round(t1/dt)+1;
i2=round(t2/dt)+1;
twin=.3;
f1=10;
f2=50;
wintype='gaus';
spectype='multi';
pflag=1;
method='specmatch';
Qint_true=Qnot;
Qint=qestimator(qmat(:,i1),qmat(:,i2),t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
prepfig
%%
%test qestimator
%a single simulation using a nonstationary reflection seismogram
dt=.002;
tmax=2;
Qnot=50;
fdom=30;
[r,t]=reflec(tmax,dt);
[w,tw]=wavemin(dt,fdom,.4,5);
Q=Qnot*ones(size(t));
flag=1;
qmat=qmatrix(Q,t,w,tw,flag);
s=qmat*r;
figure
subplot(2,1,1)
plot(t,r,t,pad(w,s)+.1,t,s+.2)
title('Nonstationary synthetic trace')
legend('reflectivity','wavelet','trace')
subplot(2,1,2)
[S,f]=fftrl(s,t);
R=fftrl(r,t);
W=fftrl(pad(w,s),t);
plot(f,real(todb(R)),f,real(todb(W)),f,real(todb(S)))

t1=.5;
t2=1.0;
twin=.3;
f1=10;
f2=100;
wintype='gaus';
spectype='multi';
pflag=1;
method='specmatch';
Qint_true=Qnot;
Qint=qestimator(s,s,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
prepfig
%% 
% an ensemble simulation using reflectivity
ntrials=200;
dt=.002;
tmax=2;
Qnot=30;
fdom=30;
[w,tw]=wavemin(dt,fdom,.4,5);
Qint=zeros(1,ntrials);
t1=.5;
t2=1.0;
twin=.5;
f1=10;
f2=50;
wintype='gaus';
spectype='fourier';
pflag=0;
method='specmatch';
% method='specrat';
Qint_true=Qnot;
t=0:dt:tmax;
Q=Qnot*ones(size(t));
driftflag=1;
qmat=qmatrix(Q,t,w,tw,driftflag);
seed=randi(1000,1,ntrials);
R=zeros(length(t),ntrials);
for k=1:ntrials
    [r,t]=reflec(tmax,dt,.1,3,seed(k));
    R(:,k)=r;
    s=qmat*r;
    Qint(k)=qestimator(s,s,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
    if(rem(k,10)==0)
        disp(['Finished trial ' int2str(k) ' of ' int2str(ntrials)])
    end
end
figure
hist(Qint,20);
title(['Q estimate, method=' method ' spectype=' spectype ', '...
    int2str(ntrials) ' trials, true value=' int2str(Qint_true)]) 
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
h=line([Qint_true Qint_true],[0 .9*ylim(2)],'color','r','linewidth',2);

seed=randi(1000,1,ntrials);
spectype='burg';    
for k=1:ntrials
    [r,t]=reflec(tmax,dt,.1,3,seed(k));
    s=qmat*r;
    Qint(k)=qestimator(s,s,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
    if(rem(k,10)==0)
        disp(['Finished trial ' int2str(k) ' of ' int2str(ntrials)])
    end
end
figure
hist(Qint,20);
title(['Q estimate, method=' method ' spectype=' spectype ', ' ...
    int2str(ntrials) ' trials, true value=' int2str(Qint_true)])
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
h=line([Qint_true Qint_true],[0 .9*ylim(2)],'color','r','linewidth',2);

spectype='multi';    
for k=1:ntrials
    [r,t]=reflec(tmax,dt,.1,3,seed(k));
    s=qmat*r;
    Qint(k)=qestimator(s,s,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true);
    if(Qint(k)==250)
        disp('hey');
    end
    if(rem(k,10)==0)
        disp(['Finished trial ' int2str(k) ' of ' int2str(ntrials)])
    end
end
figure
hist(Qint,20);
title(['Q estimate, method=' method ' spectype=' spectype ', '...
    int2str(ntrials) ' trials, true value=' int2str(Qint_true)])
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
h=line([Qint_true Qint_true],[0 .9*ylim(2)],'color','r','linewidth',2);