function phiout=ipsf(phiin,f,dx,parms,dz)
%phiout=ipsf(phiin,f,dx,parms,dz)
%
%Isotropic focussing phase shift extraploation (stationary). The this lens
%term (bulk time delay) is not applied
%
%Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005,
%Planned seismic imaging using explicit one-way operators, Geophysics, V70,
%NO5, S101 - S109' for details.
%
%phiout...f-kx spectrum* of extrapolated wavefield.
%phiin...f-kx spectrum* of input wavefield.
%f...frequency axis (Hz).
%dx...trace spacing (m).
%parms...velocity (m/s).
%dz...depth distance through which to extrapolate.
%
%*Please see calling function for details of 'spectrum'.
%
%G. F. Margrave 2014
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

%***get sizes of things***
[rp, cp]=size(phiin);
f=f(:);
rf=size(f,1);
parms=parms(:);
[rparms, cparms]=size(parms);
%*************************

%***check input***
if rparms~=1;error('  to many parmameters for isotropy');end
if cparms~=1;error('  not a stationary velocity');end
if rf~=rp;error('  frequency axis incorrect');end
%*****************

%***initialize some variables***
kx=fftshift(1/2/cp/dx*(-cp:2:cp-2));%wavenumbers (uncentered).
%*******************************

%***extrapolate one dz***
% kz=sqrt((f*ones(1,cp)/parms).^2-(ones(rp,1)*kx).^2);%vertical slowness
ff=f*ones(1,cp);
kk=ones(rp,1)*kx;
kz=(ff./parms).*(sqrt(1-(kk.*parms./ff).^2)-1);%focussing part of vertical slowness
kz=real(kz)+sign(dz)*1i*abs(imag(kz));%ensures evanescent region is complex positive
gazx=exp(2*pi*1i*dz*kz);%evenescent region will decay exponentially
phiout=phiin.*gazx;%phase shift the input spectrum
%************************
