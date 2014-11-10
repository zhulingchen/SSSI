clear
%Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005,
%Planned seismic imaging using explicit one-way operators, Geophysics, V70,
%NO5, S101 - S109' for details.

dip=70;%maximum dip (degrees).
fmax=38;%maximum frequency (Hz) of interest.
fmin=4;%minumum frequency (fmin > 0 Hz) of interest.
load Pro_salt_vel%EAGE/SEG salt model.
load Pro_salt_zero%Synthetic data (exploding reflector).

%***get sizes of things***
[rd,cd]=size(data);%data is rd rows by cd columns.
%*************************

%***build freqeuncy axis of interest***
f=1/2/dt/rd*[0:2:rd-2];%positive frequency axis in Hz.
nf=min(find(f>=fmax));%index in f for fmax.
f1=max(find(f<=fmin));%index in f for fmin.
f=f(f1:nf);%desired range of positive frequencies.
%**************************************

%***2D FFT***
temp=ifft(fft(data),[],2);clear data%fft t->f and ifft x-> kx.
fdata=temp(f1:nf,:);%spectrum band-limited between fmin and fmax.
%   NOTE f axis is band-limited and +ve only.
%   NOTE kx axis is uncentred.
%************

%***migrate***
disp('Your output file will be "fx_salt_mig"')
fx_zero=fx_zero_mig(fdata,f,model/2,dx,dz,dip);
%*************

%***output***
save fx_salt_zero_mig fx_zero
%************
