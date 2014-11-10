function mdata=ss_3c_mig(fdata,f,parms,s_parms,dx,dz)
%mdata=ss_3c_mig(fdata,f,parms,s_parms,dx,dz)
%
%***prestack splitestep 3c depth migration***
%
%To see an example, type 'ss_salt_3c_psdm_script' in matlab.
%
%Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005,
%Planned seismic imaging using explicit one-way operators, Geophysics, V70,
%NO5, S101 - S109' for details.
%
%mdata...3c depth migrated output.
%fdata...f-kx spectrum of zero-offset 'fdata=ifft(fft(data),[],2);'.
%        NOTE f axis is band-limited and +ve only.
%        NOTE kx axis is uncentred.
%        NOTE, relative to mxn velocity model 'parms', it is assumed that
%        the source is located at n/2+1.
%f...frequency axis (Hz, must correspond to fdata).
%parms...p-velocity model (m/s).
%s_parms...s-velocity model (m/s).
%dx...trace spacing (m).
%dz...depth interval (m).
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
[Nzs,Nxs]=size(s_parms);
[rd,cd]=size(fdata);
[rf,cf]=size(f(:));
%*************************

%***Err check***
if rf~=rd; error('freq axis wrong');end
if or(Nx~=cd,Nxs~=cd); error('x axis wrong');end
if or(Nxs~=Nx,Nzs~=Nz); error('p and s model dimensions must match');end
%***************

%***build the source***
temp=zeros(rd,cd);
temp(:,round(cd/2)+1)=1;
fsou=ifft(temp,[],2);
%**********************

mdata=zeros(Nz,cd);
ftemp=zeros(rd,cd);
stemp=zeros(rd,cd);
rtemp=zeros(rd,cd);
for j=1:Nz-1
	disp([' ss prestack 3c mig working on depth ',num2str(j),' of ',num2str(Nz)])
	ftemp=ss_ips(fdata,f,dx,s_parms(j,:),dz);
	fdata=fft(ftemp,[],2);
	stemp=ss_ips(fsou,f,dx,parms(j,:),-dz);
	fsou=fft(stemp,[],2);
	rtemp=fdata.*conj(fsou);%trivial reflectivity estimate
	mdata(j+1,:)=real(sum(rtemp)+sum(rtemp(1:rd-1,:)))/(2*rd-1)/2/pi;
	fdata=ftemp;
	fsou=stemp;
end
