function snew=modifysonic(s,z,t,to,tw,ts)
% MODIFYSONIC is a program that will add or subtract a bulk sonic shift 
%    from the sonic log in order to match top picks.  It is the algoritim 
%    used in the GUI program STRETCHWELL   
%
% snew=modifysonic(s,z,t,to,tw,ts)
%
%   s = sonic
%   z =depth
%   t =time vector
%   to=time at top of section where seismic and well agree
%   tw=bottom pick in well
%   ts= bottom pick in seismic
% Variables out
%   snew = New sonic log
%
% H.J.E. Lloyd, December 2012
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
dzz=z(2)-z(1);
tnew=(1/10^6)*cumsum(s.*dzz)*2;
tnew=tnew-tnew(1);
znew=interp1(tnew,z,t);
% [tz,zt,vins]=sonic2tz(s,z,-10000);
% tnew=interp1(zt,tz,z);
% znew=interp1(tz,zt,t);
alim=t(2)-t(2)/2;

io=near(t,to);
is=near(t,ts);
iw=near(t,tw);
zo=znew(io);zo=zo(1);
zs=znew(is);zs=zs(end);
zw=znew(iw);zw=zw(end);




iz1=near(z,zo);iz1=iz1(1);
iz2=near(z,zw);iz2=(iz2(end));

%%
dt=(ts-tw)/2;
a=(10^6*dt)/((zw-zo)); % assume ds=a

snew=s;
snew(iz1:iz2)=s(iz1:iz2)+a;

tt=(1/10^6)*cumsum(snew.*dzz)*2;
tt=tt-tt(1);
k=1;
tn(k)=tt(iz2);
disp(['New Time:',num2str(tn(end)),' , Number of Iterations:',num2str(k),' , Alpha:',num2str(a)]);



