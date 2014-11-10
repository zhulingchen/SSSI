function [vr,xr]=kvdesample(v,x,dxnew)
% KVDESAMPLE: desample a velocity function in the wavenumber domain
%
% vr=kvdesample(v,x,dxnew)
%
% Reduce the number of samples in a regularly sampled space series
% by reducing the maximum wavenumber. This achieves a spatial desampling
% such that the new spatial sample interval, dxr, is related to the
% original, dx, by dxr=(n/m)*dx where n is the original number of points
% and m is the new number of points. An antialias filter is applied.
%
% v ... input velocity matrix. Space series are row vectors. Each row is
%           a different depth
% x ... space coordinate vector for v
% dxnew ... new spatial sample size
% vr ... resampled velocity matrix. Will have the same number of rows as
%           v but fewer columns.

[nz,nx]=size(v);
dx=x(2)-x(1);

%determine old and new nyquist
kx=freqfft(x,[],1);
dkx=kx(2)-kx(1);
n=length(x);
knyqnew=.5/dxnew;
%knyqnew=dkx*(ceil(knyqnew/dkx)+1); % get the next largest wavenumber
ind=find(abs(kx)>knyqnew);

kxnew=kx;
kxnew(ind)=[];

m=length(kxnew);
dxnew=(n/m)*dx;

vr=zeros(nz,m);
xr=(0:m-1)*dxnew+x(1);

for k=1:nz
    tmp=fft(v(k,:));
    tmp(ind)=[];
    vr(k,:)=m*ifft(tmp)/n;
end

xr=(0:m-1)*dxnew+x(1);
