function mdata=pspi_zero_mig(fdata,f,parms,pspi_parms,dx,dz)
%mdata=pspi_zero_mig(fdata,f,parms,pspi_parms,dx,dz)
%
%***pspi zero-offset migration***
%
%To see an example, type 'pspi_salt_zero_script' in matlab.
%
%mdata...depth migrated output.
%fdata...f-kx spectrum of zero-offset 'fdata=ifft(fft(data),[],2);'.
%        NOTE f axis is band-limited and +ve only.
%        NOTE kx axis is uncentred.
%f...frequency axis (Hz, must correspond to fdata).
%parms...velocity model*1/2 (m/s).
%pspi_parms...blocky version of parms.
%dx...trace spacing (m).
%dz...depth interval (m).

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
dtemp=zeros(rd,cd);
for j=1:Nz-1
	disp([' pspi zero mig working on depth ',num2str(j),' of ',num2str(Nz)])
	dtemp=pspi_ips(fdata,f,dx,parms(j,:),pspi_parms(j,:),dz);
	mdata(j+1,:)=real(sum(dtemp)+sum(dtemp(1:rd-1,:)))/(2*rd-1)/2/pi;
	fdata=ifft(dtemp,[],2);
end
