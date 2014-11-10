function [vp,x,z,rho,I]=sandwedge(dx,xmax,zmax)
% This function creates a model with a gradient overburden and a sand wedge
% in a high velocity medium
% [vp,x,z,rho,I]=sandwedge(dx,xmax,zmax)
%
% dx ... spatial grid size (same in x and z) in physical units
% xmax ... length of model in physical units
% ******* default 3000 *******
% zmax ... depth of model in physical units
% ******* default 1000 *******
% vp ... velocity model matrix
% x ... x coordinate vector for vp
% z ... z coordinate vector for vp
% rho ... density model matrix
% I ... impedance model matrix
%
% plotting example:
% figure;
% imagesc(x,z,vp);colorbar
%
% by Heather Lloyd
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

if(nargin<3)
    zmax=1000;
end
if(nargin<2)
    xmax=3000;
end
z=0:dx:zmax;
x=0:dx:xmax;
obd=300;
chanshft=100;
chscl=.5;
% Velocity Model :
vp=3800*ones((round(zmax/dx)+1),(round(xmax/dx)+1));
vp(1:near(z,obd),:)=linspace(2000,3800,near(z,obd))'*ones(1,length(x));
vp=afd_vmodel(dx,vp,3500,[-38.313802      3022.6367      3034.3945     -38.313802     -38.313802],chanshft+chscl*[1472.9412      791.68067      1399.7479       1602.437      1472.9412]);

%Density Model :
rho=2600*ones((round(zmax/dx)+1),(round(xmax/dx)+1));
rho(1:near(z,obd),:)=linspace(1250,2600,near(z,obd))'*ones(1,length(x));
rho=afd_vmodel(dx,rho,2400,[-38.313802      3022.6367      3034.3945     -38.313802     -38.313802],chanshft+chscl*[1472.9412      791.68067      1399.7479       1602.437      1472.9412]);

I=rho.*vp;