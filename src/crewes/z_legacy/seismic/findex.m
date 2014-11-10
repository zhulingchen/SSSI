function iex=findex(trin,flag)
% iex=findex(trin,flag)
% iex=findex(trin)
%
% FINDEX returns a vector of indicies of the samples of trin 
% which are either local maxima or local minima.
% Note that extrema which persist for more than one sample will
% have only the final sample flagged.
%
% trin= input trace
% flag=1.0 .... find local maxima
%     =0.0 .... find both maxima and minima
%      -1.0 ... find local minima
%  ******* default=1.0 *******
%
% iex= vector of indicies of the extrema
%
% by G.F. Margrave June, 1991
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
if nargin<2, flag=1.0; end
 d1=diff(trin);
 ind=find(d1~=0);
 d1=d1(ind)./abs(d1(ind));
 %d1 is now +1 for a pos difference -1 for neg and zero otherwise
 d2=diff(d1);
%
% d2=-2.0 is a maximum and d2=+2.0 is a minimum
%
 
 if(flag>0.0),iex=find(d2<-1.9);end
 if(flag<0.0),iex=find(d2>1.9);end
 if(flag==0.0),iex=find(d2~=0.0);end
 iex=iex+1;
 iex=ind(iex);
