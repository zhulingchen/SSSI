function [spec,kx]=reduce_kx(inspec,v,inkx,w,dip);
%[spec,kx]=reduce_kx(inspec,vel,inkx,f,dip);
%
%Reduces the number of kx values required to extrapolate a wavefield.
%New nyquist calculated using highest frequency, slowest velocity and highest
%expected dip.
%
%  spec...spectrum of kx reduced input wavefield
%  kx...vector of reduced spatial frequencies
%  inspec...spectrum of input wavefield
%  vel...slowest scalar velocity (one way) - make sure units correspond
%        with inkx and w
%  inkx...vector of input spatial frequency 
%  w...highest temporal frequency
%  dip...highest expected dip (deg)
v=real(v);
%***compute new nyquist***
kxn=w*sin(pi*dip/180)/v;%new nyquist
dkx=abs(inkx(2)-inkx(1));
kxn=dkx*(1+fix(kxn/dkx));%ensures existance of kx=0
%*************************
%***compute new kx***
kx=[-kxn:dkx:kxn-dkx];
if abs(kx(1))>abs(inkx(1))
  clear kx;
  kx=inkx;
end
%********************
%***make a reduced spectrum***
kxs=1+(length(inkx)-length(kx))/2;
spec=inspec(:,kxs:kxs+length(kx)-1);
%*****************************
%***taper the edges***
spec(:,1)=.5*spec(:,1);
spec(:,length(kx))=.5*spec(:,length(kx));
%*********************
