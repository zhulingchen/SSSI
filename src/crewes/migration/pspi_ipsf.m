function phiout=pspi_ipsf(phiin,f,dx,parms,pspi_parms,dz)
%phiout=pspi_ipsf(phiin,f,dx,parms,pspi_parms,dz)
%
%Isotropic pspi extrapolation, focussing term only. For use in time
%migration. This is used by pspi_stack_tmig.
%
%phiout...f-x spectrum* of extrapolated wavefield.
%phiin...f-kx spectrum* of input wavefield.
%f...frequency axis (Hz).
%dx...trace spacing (m).
%parms...velocity (m/s).
%pspi_parms...low-pass filtered** velocity (m/s).
%dz...depth distance through which to extrapolate.
%stride...%1/2/dx/stride is the spatial Nyquist of the low-pass filter.
%
%*Please see function 'pspi_zero_mig.m' for details of 'spectrum'.
%
%G.F. Margrave 2014 (Modified from pspi_ips)
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
[rf, cf]=size(f);
parms=parms(:);
pspi_parms=pspi_parms(:);
% [rparms, cparms]=size(parms);
%*************************

%***check input***
if rf~=rp;error('  frequency axis incorrect');end
%*****************

%***initialize some variables***
kx=ones(rp,1)*fftshift(1/2/cp/dx*[-cp:2:cp-2]);%wavenumbers (uncentered).
k=f(:)*(1./parms');
k0=f(:)*(1./pspi_parms');
%*******************************

%***extrapolate one dz***
vref=unique_vels(pspi_parms);
rv=length(vref(:));
phiout=zeros(rp,cp);
for j=1:rv
	inds=find(pspi_parms==vref(j));
	temp1=fft(ipsf(phiin,f,dx,vref(j),dz),[],2);
	phiout(:,inds)=phiout(:,inds)+temp1(:,inds);
end
% thin_lens=exp(2*pi*1i*dz*(k-k0));
% phiout=thin_lens.*phiout;
%************************
