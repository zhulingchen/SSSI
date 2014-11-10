function [x,t,phi]=rayfan(v,z,z1,z2,p)
%
% [x,t,phi]=rayfan(v,z,z1,z2,p)
%
% rayfan shots a fan of rays through the stratified earth. It is
% identical to rayfan_a except that the fan is described by a p vector
% rather than an angle vector.
% v = vector of layer velocities
% z = vector of depths to the tops of each layer. Thus V(k) is
%	the velocity between z(k) and z(k+1). Z should be sorted into
%	ascending order.
% z1 = depth the ray starts at
% z2 = depth the ray ends at
% p = vector of ray parameters defining the fan
% x = horizontal distance the ray travels (one for each p)
% t = traveltime along the ray ( in seconds ) (one for each p)
% NOTE: Critically refracted rays will return inf for both x an t
% phi = angle (from the vertical) of the ray at z2 (degrees) (one for each p)
%
% G.F. Margrave, CREWES Project, July 1995
%preliminaries
 z=z(:);
 v=v(:);
 p=p(:)'; %make p a row vector
 nz=length(z);
 %adjust z1 and z2 so that they won't be exactly on layer boundaries
 zt=min([z1 z2]);
 if(zt==z1)
 	z1=z1+100000*eps;
 	z2=z2-100000*eps;
 else
 	z1=z1-100000*eps;
 	z2=z2+100000*eps;
 end
 if( z1< z(1) | z2< z(1) )
 	error(' start or end depth outside model range');
 end
%determine layers propagated through
ind=find(z>=z1);
if(isempty(ind))
 	ibeg=nz;
 else
 	ibeg=ind(1)-1;
 end
 
 %compute sin and check for critical angle
 ind=find(z>=z2);
 if(isempty(ind))
 	iend=nz;
 else
 	iend=ind(1)-1;
 end
 if(ibeg~=iend)
 	 % these are the layers propagated thru
 	iprop=ibeg:sign(iend-ibeg):iend;
 else
 	iprop=ibeg;
 end
 sn = v(iprop)*p;
 %
 % sn is an n by m matrix where n is the length of iprop (the number of
 %    layers propagated through) and m is the length of p (the number of
 %    unique ray parameters to use). Each column of sn corresponds to a
 %    single ray parameter and contains the sin of the vertical angle;
 %
 
 [ichk,pchk]=find(sn>1);
% if(~isempty(ichk))	
% 	%assign zeros for critical refractions
% 	sn(ichk,pchk)=zeros(size(sn(ichk,pchk)));
% end
 
 %compute x and t
 np=length(iprop);
 cs=sqrt(1-sn.*sn);
 zprop=[z1;z(iprop(2:np));z2];
 vprop=v(iprop)*ones(1,length(p));
 thk=abs(diff(zprop))*ones(1,length(p));
 if(size(sn,1)>1)
 	x=sum( (thk.*sn)./cs);
 	t=sum(thk./(vprop.*cs));
 else
 	x=(thk.*sn)./cs;
 	t=thk./(vprop.*cs);
 end
 if(~isempty(ichk))
 	x(pchk)=inf*ones(size(pchk));
 	t(pchk)=inf*ones(size(pchk));
 end
 
 % final angle
 if(nargout>2)
 	phi=180*asin(sn(np,:))/pi;
 end
