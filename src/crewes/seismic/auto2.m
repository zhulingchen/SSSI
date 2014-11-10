function a = auto2(v,flag)
% AUTO2: returns the two-sided autocorrelation 
%
% a=auto2(v,flag)
% a=auto2(v)
%
% auto computes the two sided autocorrelation of 
% the vector 'v'. The 'zeroth lag' is at length(v)-1. This 
% routine will correctly handle an ensemble matrix.
% 
% v= input vector
% flag= 1.0 ... normalize the 'zero lag' (first returned value)
%               to 1.0.
%        anything else ... don't normalize
%       ******* default =1.0 ******
%
% NOTE: A two sided autocorrelation or a cross correlation can be
%       computed with XCORR in the signal toolbox.
%
%   by G.F. Margrave, July 1991
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
 
% set defaults
 if (nargin<2)
  flag=1.0;
 end
% 
% 
 nvecs=min(size(v));
 if nvecs==1,
  a=conv(v,tr(v,1));
% normalize
  if flag==1.0
    a=a/max(a);
  end
 else
  [rows,cols]=size(v);
  if rows>cols, v=v.'; [rows,cols]=size(v); end
  a=zeros(rows,2*cols-1);
  if flag==1,
   for k=1:rows,
     a(k,:)=conv(v(k,:),tr(v(k,:),1));
     a(k,:)=a(k,:)/max(a(k,:));
   end
  else
   for k=1:rows,
     a(k,:)=conv(v(k,:),tr(v(k,:),1));
   end
  end
 end
			     
