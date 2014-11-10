function imp=rcs2imp_r(rcs,znot,flag)
% imp=rcs2imp_r(rcs,znot,flag)
% imp=rcs2imp_r(rcs,znot)
%
% RCS2IMP_R computes impedence from RCS by recursion. It is
% more accurate but slower than RCS2IMP
% 
% rcs ... input column vector of reflection coefficients
% znot ... scalar value of first impedence
% flag ... if 1, then a shift of 1/2 a sample interval is
%	applied. Otherwise, no action. This is useful to more
%	properly invert the rc computation in imp2rcs
% ******** default = 0 ***************
%
% G.F. Margrave May 1995
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
if(nargin<3) flag=0; end
n=length(rcs);
imp=zeros(size(rcs));
imp(1)=znot;
for k=2:n
	imp(k)=imp(k-1)*(1+rcs(k-1))/(1-rcs(k-1));
end
if(flag==1)
	%make a time coordinate vector(arbitrary)
	t=xcoord(0.,.001,rcs);
	%fft
	[Imp,f]=fftrl(imp,t);
	%design phase shifter
	phs=cos(pi*.001*f)+i*sin(pi*.001*f);
	%apply
	Imp=Imp.*phs;
	%inverse fft
	[imp,t]=ifftrl(Imp,f);
end
