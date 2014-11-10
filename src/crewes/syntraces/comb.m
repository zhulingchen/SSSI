function trout=comb(trin,n,flag)
%COMB ... create a comb function (spikes every n samples)
%
%    trout=comb(trin,n,flag) 
%    trout=comb(trin,n)
%
% outputs a comb function which has unit spikes every n samples
% 
% trin= input trace (to give dimensions)
% n= interval between spikes in samples
%    if n is negative, the spikes will alternatly be +1 and -1
% flag= 0 ... first unit spike is at sample n/2 (*** default ***)
% flag= 1 ... first unit spike is at sample 1
% flag= 2 ... first unit spike is at sample n
% trout= output comb function
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

if nargin<3
   flag=0;
 end
%
   trout=zeros(size(trin));
   m=floor(abs(n/2));
   if flag==1
     m=1;
   end
   if flag==2
     m=abs(n);
   end
   if n>0
     while m<=length(trin)
       trout(m)=1.0;
       m=m+n;
     end
   else
     spike=-1.0;
     while m<=length(trin)
       spike=-1*spike;
       trout(m)=spike;
       m=m+abs(n);
     end
   end
