%make a section full of diffractions
dt=.002;t=((0:500)*dt)';
dx=10;x=(0:200)*dx;
fmin=[10 5];fmax=[125 25];
nt=length(t);nx=length(x);
seis=zeros(nt,nx);
%determine diffraction locations
itdiff=1:50:nt;
ixdiff=1:20:nx;    
tmp=zeros(size(t));
tmp(itdiff)=1;
%bandlimit
tmp=filtf(tmp,t,fmin,fmax);
%install
for k=ixdiff
    seis(:,k)=tmp;
end
%linear v(z) function
z=0:dx:2000;
vz=1800+.6*z;
tv=2*vint2t(vz,z);
v=interp1(tv,vz,t);
%constant velocity
vm=mean(v);
vc=ones(size(v))*vm;
%set up for modelling
params=nan*(1:13);
params(1)=150;
params(2)=0;
params(13)=1;
[seismod,tmod,xmod]=vz_fkmod(seis,v,t,x,params);
plotimage(seismod,tmod,xmod);
[seismodc,tmod,xmod]=vz_fkmod(seis,vc,t,x,params);
plotimage(seismodc,tmod,xmod);