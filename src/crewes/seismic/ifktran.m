function [seis,t,x]=ifktran(spec,f,kx,nfpad,nkpad,percent)
% IFKTRAN inverse fk transform
% [seis,t,x]=ifktran(spec,f,kx,nfpad,nkpad,percent)
%
% IFKTRAN uses Matlab's built in fft to perform a 2-d inverse f-k transform
% on a complex valued (presumably seismic f-k) matrix.  The forward transform
% is performed by fktran.
%
% spec ... complex valued f-k transform
% f ... vector of frequency coordinates for the rows of spec.  
%	length(f) must be the same as number of rows in spec.
% kx ... vector of wavenumber coordinates for the columns of spec
%	length(kx) must be the same as the number of columns in spec.
% nfpad ... Non-functional, enter 0 to specify later parameters.
%	******* default = 0 ******
% nkpad ... pad spec with zero filled columns until it is this size.
%   ******* default = 0 ******
% NOTE: a value of 0 for nkpad is taken to mean not padding
% percent ... apply a raised cosine taper to both f and kx prior to zero pad.
%	length of taper is theis percentage of the length of the kx and f axes
%   ******* default = 0 *********
% seis ... output 2-d seismic matrix. One trace per column.
% t ... vector of time coordinates for seis. 
% x ... vector of space coordinates for seis. 
% 
% G.F. Margrave, CREWES Project, U of Calgary, 1996
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

[nf,nkx]=size(spec);

if(length(f)~=nf)
	error(' Frequency coordinate vector is incorrect');
end
if(length(kx)~=nkx)
	error(' Wavenumber coordinate vector is incorrect');
end

if(nargin<6) percent=0.; end
if(nargin<5) nkpad=0; end
if(nargin<4) nfpad=0; end

%determine if kx needs to be wrapped
if(kx(1)<0) % looks unwrapped
	ind=find(kx>=0);
	kx=kx([ind 1:ind(1)-1]);
	spec=spec(:,[ind 1:ind(1)-1]);
end

% ok taper and pad in kx
if(percent>0)
	mw=mwindow(nkx,percent);
	mw=mw([ind 1:ind(1)-1]);
	mw=mw(ones(1,nkx),:);
	spec=spec.*mw;
	clear mw;
end

if(nkx<nkpad)
	nkx=nkpad;
end

%ok kx-x transform
%disp('kx-x')
specfx= fft(spec.',nkx).';
clear spec;

%use ifftrl for the f-t transform
%disp('f-t');
[seis,t]=ifftrl(specfx,f);
clear specfx;

% compute x
dkx=kx(2)-kx(1);
xmax = 1./(dkx);
dx = xmax/nkx;
x=0:dx:xmax-dx;


