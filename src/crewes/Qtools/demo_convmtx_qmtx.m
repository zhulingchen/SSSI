% make a synthetic reflectivity and a minimum phase wavelet
dt=.002;
tmax=2;
fdom=30;
[r,t]=reflec(tmax,dt);%synthetic reflectivity
[w,tw]=wavemin(dt,fdom,.2);%min phase wavelet
s=convm(r,w);%convolve them using convm
%now plot the results using subplot to make 3 rows in the figure
figure
subplot(3,1,1)
plot(tw,w);%note this wavelet plot will have a different time scale from the other subplots
subplot(3,1,2)
plot(t,r)
subplot(3,1,3)
plot(t,s)
%now lets fix up the first subplot to have the same time scale as the
%others. We us "pad" for this
subplot(3,1,1)
plot(t,pad(w,r))
% now lets do the convolution with a convolution matrix
%build a convolution matrix from w that is the right size to convolve with
%r
cmtx=convmtx(w,length(r));%the convolution matrix
plotimage(cmtx)%view it with plotimage
title('Convolution matrix with all rows')
cmtxm=cmtx(1:1001,:);%grabbing the first 1001 rows because we want to simulate convm
plotimage(cmtxm)
title('Convolution matrix after row truncation')
s2=cmtxm*r;%this should be identical to s. Is it?
figure
plot(t,s,t,s2,'r.')%this plot suggests s and s2 are pretty close
legend('Convm','Convolution matrix')
a=sum(abs(s-s2))/length(s);%this is a more precise test of equivalence
b=eps;%eps tells me thee precision of my computer
title(['Convm compared to convmtx, ave error=' num2str(a) ' and eps=' num2str(b)])
% the fact that sum(abs(s-s2)) is similar to eps tells me that s and s2
% are equivalent to machine precision.
%%
%another way to make a convolution matrix is with the qmatrix command
%Q is the measure of attenuation in rocks. A Q of infinity means no
%attenuation while a Q of 50 is lots of attenuation.
Q=50;
qmat=qmatrix(inf,t,w,tw);%inf means "infinity" here.
plotimage(qmat)
title('Stationary Q matrix')
qmat2=qmatrix(Q,t,w,tw);%another matrix for Q=50.
plotimage(qmat2)
title(['Nonstationary Q matrix, Q=' int2str(Q)])
%compare the above two Q matrices. Note the difference between stationary
%and nonstationary
sinf=qmat*r;%stationary trace
sQ=qmat2*r;%nonstationary trace
figure
plot(t,sinf,t,sQ,'r')%compare the traces in the time domain
legend('Stationary',['Nonstationary, Q=' int2str(Q)])
prepfig
xlabel('Time (s)')
figure
hh=dbspec(t,[sinf sQ])%compare the traces in the frequency domain;
set(hh(2),'color','r')
legend('Stationary',['Nonstationary, Q=' int2str(Q)])
prepfig
%%
%examine local spectra in two windows
t1=.2;
t2=1.2;
twin=.6;
inwin1=near(t,t1,t1+twin);%finds the samples in window 1
inwin2=near(t,t2,t2+twin);%finds the samples in window 2
[S1inf,f]=fftrl(sinf(inwin1),t(inwin2));%Spectrum in win1
S2inf=fftrl(sinf(inwin2),t(inwin2));%spectrum in win2
S1Q=fftrl(sQ(inwin1),t(inwin2));%Spectrum in win1
S2Q=fftrl(sQ(inwin2),t(inwin2));%spectrum in win2
figure
plot(f,abs(S1inf),f,abs(S2inf),'r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Stationary case')
legend(['Window 1 at ' num2str(t1+.5*twin) ' s'],...
    ['Window 2 at ' num2str(t2+.5*twin) ' s'])
prepfig
figure
plot(f,abs(S1Q),f,abs(S2Q),'r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title(['Nonstationary case, Q=' int2str(Q)])
legend(['Window 1 at ' num2str(t1+.5*twin) ' s'],...
    ['Window 2 at ' num2str(t2+.5*twin) ' s'])
prepfig
%%
[Qmat2,f]=fftrl(qmat2,t);%nonstationary
[Qmat,f]=fftrl(qmat,t);%stationary
figure
subplot(1,2,1)
imagesc(t,f,abs(Qmat))
ylim([0 100])
xlabel('time (sec)')
ylabel('frequency (Hz)')
title('Stationary')
subplot(1,2,2)
imagesc(t,f,abs(Qmat2))
ylim([0 100])
xlabel('time (sec)')
title('Nonstationary')
ylabel('frequency (Hz)')
prepfig