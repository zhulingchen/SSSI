function L=near(v,val1,val2)
% L=near(v,val1,val2)
% L=near(v,val1)
%
% NEAR searches the vector v and finds the index,L1, for which
%   v(l) is closest to val1 and L2 for which v(i) is closest to
% val2. The returned valued is the vector L=L1:L2 (L2>L1) or
% L=L1:-1:L2 for L2<L1 
%
% v= input vector
% val1= first input search value
% val2= second input serach value
%  ******** default= val1 ********
% L= output vector of indicies such that 
% abs(v(l(1))-val1)==minimum and abs(v(l(length(I))-val2)==minimum
%
% by G.F. Margrave, May 1991
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
  
% 
 if nargin<3
    val2=val1;
 end
 
  ilive=find(~isnan(v));
  test=abs(v(ilive)-val1);
  L1=find(test==min(test));
  test=abs(v-val2);
  L2=find(test==min(test));
  L1=ilive(L1);
  L2=ilive(L2);
 	if L1<=L2
  		L=min(L1):max(L2);
	else
    	L=max(L1):-1:min(L2);
 	end
 
