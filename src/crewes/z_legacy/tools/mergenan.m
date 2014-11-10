function vout=mergenan(v1,v2,w1,w2)
% vout=mergenan(v1,v2,w1,w2)
%
% MERGENAN forms the weighted sum of 2 matricies with special attention 
% paid to nan's. If both inputs are nan at any position, then the output
% will also be; however, if only one is nan then the output will equal
% the live matrix at that location. If both are live then the output is
% the direct weighted sum. If v1 and v2 contain no nans, then the result 
% will be exactly ((w1.*v1)+(w2.*v2))./(w1+w2); 
% The weights can be simple scalars or they can be matricies of exactly
% the same size as v1 and v2.
%
% G.F. Margrave March 1994
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
i1nan=isnan(v1);
i2nan=isnan(v2);
[m,n]=size(w1);
if( m==1 & n==1)
	w1=w1*ones(size(v1));
end
[m,n]=size(w2);
if( m==1 & n==1)
	w2=w2*ones(size(v2));
end
bothlive= (~i1nan).*(~i2nan);
vout=nan*ones(size(v1));
aa=find(bothlive==1);
vout(aa)= ((w1(aa).*v1(aa))+(w2(aa).*v2(aa)))./...
(w1(aa)+w2(aa));
	
i1dead=find(i1nan);
i2live=find(~isnan(v2(i1dead)));
vout(i1dead(i2live))=v2(i1dead(i2live));
i2dead=find(i2nan);
i1live=find(~isnan(v1(i2dead)));
vout(i2dead(i1live))=v1(i2dead(i1live));
