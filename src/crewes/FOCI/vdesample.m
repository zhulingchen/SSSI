function [vr,xr]=vdesample(v,x,dxnew)
% VDESAMPLE: desample a velocity function by spatial averaging
%
% vr=vdesample(v,x,dxnew)
%
% Reduce the number of samples in a regularly sampled space series
% by spatial averaging
%
% v ... input velocity matrix. Space series are row vectors. Each row is
%           a different depth
% x ... space coordinate vector for v
% dxnew ... new spatial sample size
% vr ... resampled velocity matrix. Will have the same number of rows as
%           v but fewer columns.
%
% G.F. Margrave, CREWES/POTSI 2004
%


[nz,n]=size(v);
dx=x(2)-x(1);

%m=round((max(x)-min(x))/dxnew)+1;
m=round(n*dx/dxnew);%should be integer, round just in case
xr=(0:m-1)*dxnew+x(1);

vr=zeros(nz,m);

for k=1:m-1
    ind=between(xr(k),xr(k+1),x,1);
    vr(:,k)= mean(v(:,ind),2);
end
vr(:,m)=vr(:,m-1);
