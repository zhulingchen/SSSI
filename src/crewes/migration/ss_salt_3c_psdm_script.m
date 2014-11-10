clear
%Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005,
%Planned seismic imaging using explicit one-way operators, Geophysics, V70,
%NO5, S101 - S109' for details.

fmax=38;%maximum frequency (Hz) of interest.
fmin=4;%minumum frequency (fmin > 0 Hz) of interest.
load Pro_salt_vel%EAGE/SEG salt model.
[rm cm]=size(model);
s_model=model/2;
load shot_impulse
dt=t(2);
dx=xrec(2);
dz=4;

%***get sizes of things***
[rs,cs]=size(seis);%seis is rs rows by cs columns.
%*************************

%***build freqeuncy axis of interest***
f=1/2/dt/rs*[0:2:rs-2];%positive frequency axis in Hz.
nf=min(find(f>=fmax));%index in f for fmax.
f1=max(find(f<=fmin));%index in f for fmin.
f=f(f1:nf);%desired range of positive frequencies.
%**************************************

%***2D FFT***
temp=ifft(fft(seis),[],2);clear seis%fft t->f and ifft x-> kx.
fseis=temp(f1:nf,:);%spectrum band-limited between fmin and fmax.
%   NOTE f axis is band-limited and +ve only.
%   NOTE kx axis is uncentred.
%************

%***migrate***
disp('Your output file will be "ss_3c_salt_mig"')
pause(3)
ss_3c_psdm=ss_3c_mig(fseis,f,model,s_model,dx,dz);
%*************

%***output***
save ss_salt_3c_mig ss_3c_psdm
%************
