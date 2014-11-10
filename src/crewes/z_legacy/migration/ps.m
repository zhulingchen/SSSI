function out=pstime(inspec,vel,inkx,f,dz,dip)
%
%out=ps(inspec,vel,inkx,f,dz,dip);
%
%Advance wavefield one depth step by constant velocity phase shift.
%
%  out...spectrum of extrapolated wavefield
%  inspec...spectrum of input wavefield
%  vel...scalar velocity - make sure units correspond
%        between x, dz, inkx, f
%  inkx...spatial frequencies
%  f...temporal frequencies
%  dz...depth through which to extrapolate
%  dip...dip cut off (degrees, used to limit kx - not a dip filter)
%***kx limit the input data***
[rf cf]=size(f);
[spec,kx]=reduce_kx(inspec,vel,inkx,f(rf),dip);
%*****************************
%***get sizes of things***
[rins cins]=size(inspec);clear inspec;
[rs cs]=size(spec);
%*************************
%***initialize some variables***
vel2=vel^2;
kxkx=kx.*kx;
phiout=zeros(rs,cs);
kz=zeros(1,cs);
gazx=zeros(1,cs);
%*******************************
%***extrapolate each frequency***
for j=2:rf
  kz=sqrt((f(j)*f(j)/vel2)-kxkx);
  kz=real(kz)+i*abs(imag(kz));%Evanecent inverter
  gazx=exp(2*pi*i*dz*kz);
  phiout(j,:)=spec(j,:).*gazx;
end
clear kz; clear gazx; clear spec;
%********************************
%***restore the dip limited spectrum to full size***
kxs=1+(cins-cs)/2;
out=zeros(rins,cins);
out(:,kxs:kxs+cs-1)=phiout;
%***************************************************
