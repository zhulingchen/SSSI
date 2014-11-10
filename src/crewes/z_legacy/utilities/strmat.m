function smat=strmat(s1,s2)
% smat=strmat(s1,s2)
%
% STRMAT can be used to form string matricies. It differs from MATLAB's
% STR2MAT mainly in that it pads strings with 1's instead of blanks.
% This makes detection of the pad easier when the string itself may
% contain blanks. 1's are non-displaying ASCII characters.
% If s1 is a string matrix (each separate string is a row) with n rows,
% then smat will have n+1 rows and as many columns as the larger of
% s1 & s2.
% See also STRPAD STRUNPAD
%
% G.F. Margrave, March 1994
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
	[nrow,ncol]=size(s1);
	[n,m]=size(s2);
 if( n~=1 & m==1 ) s2=s2'; end
 % determine pad 
 if( m>ncol )
		s1=[s1 ones(nrow,m-ncol)];
	elseif( ncol>m )
		s2=[s2 ones(n,ncol-m)];
	end
	smat=[s1;s2];
