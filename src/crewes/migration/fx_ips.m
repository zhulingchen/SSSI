function phiout=fx_ips(phiin,f,dx,parms,dz,dip)
%phiout=fx_ips(phiin,f,dx,parms,dz,dip)
%
%Isotropic 80 degree f-x extraploation.
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
%dip...maximum dip (usually dip < 80 degrees)
%
%*Please see calling function 'fx_zero_mig.m' for details of 'spectrum'.
%**Please see calling function 'fx_zero_mig.m' for details of 'low_pass'.
%
%R. J. Ferguson, 2009
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
[rp cp]=size(phiin);
f=f(:);
[rf cf]=size(f);
parms=parms(:)';
[rparms cparms]=size(parms);
%*************************

%***check input***
if rparms~=1;error('  to many parmameters for isotropy');end
if cparms~=cp;error('  a stationary velocity');end
if rf~=rp;error('  frequency axis incorrect');end
%*****************

%***initialize some variables***
kx=ones(rp,1)*fftshift(1/2/cp/dx*[-cp:2:cp-2]);%wavenumbers (uncentered).
k=f*(1./parms);
kmax=f/max(parms)*sin(pi*dip/180)*ones(1,cp);
%*******************************

%***extrapolate one dz***
x=-(kx./kmax).^2;
temp=sqrt(1+x);
pass=find(imag(temp));
phiin(pass)=0;
g0=1;
g1=-i*pi*dz./k;
g2=(-1/4*i*pi*dz*k-1/2*pi^2*dz^2*k.^2)./k.^4;
g3=-(1/8*i*pi*dz*k+1/4*pi^2*dz^2*k.^2-1/6*i*pi^3*dz^3*k.^3)./k.^6;
g4=(-5/64*i*pi*dz*k-5/32*pi^2*dz^2*k.^2+1/8*i*pi^3*dz^3*k.^3+1/24*pi^4*dz^4*k.^4)./k.^8;
t0=g0.*fft(phiin,[],2);
t1=g1.*fft(phiin.*kx.^2,[],2);
t2=g2.*fft(phiin.*kx.^4,[],2);
t3=g3.*fft(phiin.*kx.^6,[],2);
t4=g4.*fft(phiin.*kx.^8,[],2);
phiout=ifft(exp(2*i*pi*dz*k).*(t0+t1+t2+t3+t4),[],2);
%************************
