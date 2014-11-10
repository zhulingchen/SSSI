function poleplot(x,y,z,zref)
% poleplot( x, y, z, zref)
%
% x = vector containing the x coordinates of the data
%				 ******* default= no field set *********
% y = vector containing the y coordinates of the data
%				 ******* default= no field set *********
% z = vector containing the z coordinates of the data
%				 ******* default= no field set *********
% zref = reference zlevel for drawing x,y paths
% ************* default = mean(z) *******************
%
% by G.F. Margrave, March 1993
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE
if(nargin <4)
	zref = mean(z);
end
x=x(:);
y=y(:);
z=z(:);
n=length(z);
displayx=zeros(3,n);
displayy=zeros(3,n);
displayz=zeros(3,n);
xbar= mean(x);
ybar=mean(y);
for j=1:n
  displayx(1,j)=xbar;
  displayx(2,j)=x(j);
  displayx(3,j)=x(j);
  displayy(1,j)=ybar;
  displayy(2,j)=y(j);
  displayy(3,j)=y(j);
  displayz(1,j)=zref;
  displayz(2,j)=zref;
  displayz(3,j)=z(j);
end
plot3(displayx,displayy,displayz,'-b');
hold on;
index=find(z>zref);
plot3(x(index),y(index),z(index),'m*');
index=find(z<=zref);
plot3(x(index),y(index),z(index),'r*');
hold off;
