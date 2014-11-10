function trout=phsrot(trin,theta)
% PHSROT Constant-phase rotate a trace
% trout=phsrot(trin,theta)
%
% PHSROT performs a constant phase rotation of the input trace
% through an angle of theta degrees.
%
% trin= input trace
% theta= phase rotation angle in degrees
% trout= phase rotated trace
%
% by G.F. Margrave, June 1991
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

% make fake time coordinate
 [m,n]=size(trin);
 %trin=padpow2(trin(:).');
 t=xcoord(0.,.002,trin);
 [Trin,f]=fftrl(trin,t);
 op=exp(i*theta*pi/180.);
 trout=ifftrl(Trin*op,f);

 % reshape trout
 if(m==1)
		trout=trout(1:n);
	else
		trout=trout(1:m).';
	end


