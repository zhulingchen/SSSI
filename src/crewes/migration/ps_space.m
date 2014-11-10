function out=ps_space(spec,f,dx,v,dz)
%out=ps_space(spec,f,dx,v,dz)
%
%This function extrapolates an input Fourier spectrum through dz
%by space domain phase screen.
%
%out...x-w spectrum of extrapolated seismic data
%spec...x-w spectrum of input seismic data
%f...frequency axis of spec in Hz (cycles/second)
%dx...x interval
%v...v(x) (must correspond to the range of x spanned by the seismic data)
%dz...depth step
%
%Rob Ferguson 2003
%Assistant Professor
%Department of Geological Sciences
%University of Texas, Austin
%512 471 6405
%fergusonr@mail.utexas.edu

%***find sizes of things***
f=f(:);%force f to be a column vector
v=v(:)';%force v to be a row vector
[rs cs]=size(spec);
[rf cf]=size(f);
[rv cv]=size(v);
[rz cz]=size(dz);
%***************************

%***check input***
if rs~=rf;disp('ERROR: something is wrong with the size of the input spectrum and/or f axis');end
if cs~=cv;disp('ERROR: the velocity has to have as many elements as the spectrum is wide');end
if and(rz~=1,cz~=1);disp('ERROR: dz should be a scalar');end
%*****************

%***innitialize some variables***
vbar=mean(v(:));
%********************************

%***extrapolate each monochromatic (in w) set of phases***
temp1=exp(i*2*pi*dz*(f(:)*ones(1,cv)*vbar./((ones(rf,1)*v).^2)-f(:)*ones(1,cv)/vbar)).*spec;
temp2=ifft(temp1,[],2);
spec=ips(temp2,f,dx,vbar,dz);
out=fft(spec,[],2);
%*********************************************************
