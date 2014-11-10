function x = fftshiftc(x)
% x= fftshiftc(x)
%
% This version of FFTSHIFT supplied in the Matlab toolbox is designed
% to shift the columns of a large matrix, much as FFTSHIFTM does
% but without the large memory requirements. The operation takes
% place on each column, one at a time. For matricies of normal
% size use FFTSHIFTM first, if memory problems develop, then try
% this. For vectors, use FFTSHIFT or FFTSHIFTM.
%
% Correction by G.F. Margrave, May 1991
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
[m,n] = size(x);
m1 = 1:ceil(m/2);
m2 = ceil(m/2)+1:m;
% Note: m/2+1 references the Nyquist. so the spectrum after FFTSHIFT 
% goes from -Nyquist -> DC -> one sample before + Nyquist
icol=1;
 while icol<=n
  temp=x(:,icol);
  temp(:)=[temp(m2) temp(m1)];
  x(:,icol)= temp;
  icol=icol+1;
 end  
