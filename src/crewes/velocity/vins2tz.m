function tzcurve=vins2tz(vins,z,nlegs,tnot)
% VINS2TZ: given instantaeous velocity versus depth compute a two-way time-depth curve
% 
% tzcurve=vins2tz(vins,z,nlegs,tnot)
%
% vins ... column vector in instantaneous velocities
% z ... vector of depths to pair with vins 
% nlegs ... number of points desired on the tz curve
%  ******* default nlegs = length(z) ********  
% tnot ... two-way traveltime to depth z(1)
%  ******* default tnot=0 ********
%
% tzcurve ... matrix of dimension nlegs by 2, first column is two-way
% vertical traveltimes and the second column is the corresponding depths
% The depths are evenly sampled as linspace(min(z),max(z),nlegs)
%

if(nargin<4)
    tnot=0.0;
end
if(nargin<3)
    nlegs=length(z);
end

if(nlegs<length(z))
    zout=linspace(min(z), max(z), nlegs);
else
    zout=z;
end

t=2*vint2t(vins,z,zout,tnot/2);%make sure it's two-way time

tzcurve=[t(:) zout(:)];
    