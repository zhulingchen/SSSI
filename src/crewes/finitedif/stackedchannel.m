function [x,z,vp,rho,I]=stackedchannel
% this function creates a stacked channel model with a gradient overburden.
% [x,z,vp,rho,I]=stackedchannel
%
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
dxz=2;
zmax=1000;
xmax=3000;
z=0:dxz:zmax;
x=0:dxz:xmax;
obd=600;
chanshft=600;
chscl=.15;
% Velocity Model :
vp=3800*ones((round(zmax/dxz)+1),(round(xmax/dxz)+1));
vp(1:near(z,obd),:)=linspace(2000,3800,near(z,obd))'*ones(1,length(x));
vp=afd_vmodel(dxz,vp,3500,[428.07943      1411.8164      2524.8893      2246.6211      1968.3529      1690.0846      1411.8164      1153.1445      919.94792      674.01367      428.07943],chanshft+chscl*[904.28571      932.43697      864.87395      1046.6561      1167.6589      1238.6139      1270.2521      1237.8782      1157.0246      1002.4926      904.28571]);
vp=afd_vmodel(dxz,vp,3500,[992.45443      1972.2721      3089.2643      2810.9961      2532.7279      2254.4596      1976.1914      1717.5195      1484.3229      1238.3887      992.45443],chanshft+chscl*[617.1429      628.4034      577.7311      759.5132      880.5161       951.471      983.1092      950.7353      869.8817      715.3498      617.1429]);

%Density Model :
rho=2600*ones((round(zmax/dxz)+1),(round(xmax/dxz)+1));
rho(1:near(z,obd),:)=linspace(1250,2600,near(z,obd))'*ones(1,length(x));
rho=afd_vmodel(dxz,rho,2400,[428.07943      1411.8164      2524.8893      2246.6211      1968.3529      1690.0846      1411.8164      1153.1445      919.94792      674.01367      428.07943],chanshft+chscl*[904.28571      932.43697      864.87395      1046.6561      1167.6589      1238.6139      1270.2521      1237.8782      1157.0246      1002.4926      904.28571]);
rho=afd_vmodel(dxz,rho,2400,[992.45443      1972.2721      3089.2643      2810.9961      2532.7279      2254.4596      1976.1914      1717.5195      1484.3229      1238.3887      992.45443],chanshft+chscl*[617.1429      628.4034      577.7311      759.5132      880.5161       951.471      983.1092      950.7353      869.8817      715.3498      617.1429]);

I=rho.*vp;