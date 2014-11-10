function mdata=fx_zero_mig(fdata,f,parms,dx,dz,dip)
%mdata=fx_zero_mig(fdata,f,parms,dx,dz,dip);
%
%***fx zero offset migration***
%
%To see an example, type 'fx_salt_zero_script' in matlab.
%
%Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005,
%Planned seismic imaging using explicit one-way operators, Geophysics, V70,
%NO5, S101 - S109' for details.
%
%mdata...depth migrated output.
%fdata...f-kx spectrum of zero-offset 'fdata=ifft(fft(data),[],2);'.
%        NOTE f axis is band-limited and +ve only.
%        NOTE kx axis is uncentred.
%f...frequency axis (Hz, must correspond to fdata).
%parms...velocity model*1/2 (m/s).
%dx...trace spacing (m).
%dz...depth interval (m).
%dip...maximum dip (degrees)
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
[Nz,Nx]=size(parms);
[rd,cd]=size(fdata);
[rf,cf]=size(f(:));
%*************************

%***Err check***
if rf~=rd; error('freq axis wrong');end
if Nx~=cd; error('x axis wrong');end
%***************

mdata=zeros(Nz,cd);
temp=zeros(rd,cd);
for j=1:Nz-1
	disp([' fx zero mig working on depth ',num2str(j),' of ',num2str(Nz)])
	temp=fx_ips(fdata,f,dx,parms(j,:),dz,dip);
	fdata=ips(temp,f,dx,max(parms(j,:))/sin(pi*dip/180),dz);%stabilize
	temp=ips(fdata,f,dx,max(parms(j,:))/sin(pi*dip/180),-dz);%stabilize
	fdata=fft(temp,[],2);
	mdata(j+1,:)=real(sum(fdata)+sum(fdata(1:rd-1,:)))/(2*rd-1)/2/pi;
	fdata=temp;
end
