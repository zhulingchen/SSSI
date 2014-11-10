function [mfilt,tm]=match(trin,trdsign,t,mlength,flag)
% [mfilt,tm]=match(trin,trdsign,t,mlength,flag)
%
% MATCH designs a match filter of temporal length 'mlength'
% which matches trin to trdsign in the least squares sense.
% That is sum_of_sqs(conv(mfilt,trin)-trdsign)==minimum
%
% trin= input trace to be matched to trdsign
% trdsign= input trace which is to be matched
% t= time coordinate vector for trin
% ***** note: trin and trdsign must be the same length
% mlength= length of the match filter in seconds
% flag=0 ... a noncausal operator is desired
%     =1 ... a causal operator is desired
% NOTE: Suppose that 'a' is a known time series and 'w' is also
%  a known wavelet. Then let bm=convm(a,w) and bz=convz(a,w).
%  Then [westm,tw]=match(a,bm,ta,length(w),1) or
%       [westz,tw]=match(a,bz,ta,length(w),0) 
%  will both produce good estimates of w. 
%  However, 
%       [westx,tw]=match(a,bz,ta,length(w),1) 
%  should not be used
%  but
%       [westy,tw]=match(a,bm,ta,2*length(w),0)
%  will produce a valid estimate of w in the second half of westy.         
%
% mfilt= output mlength match filter
% tm= time coordinate for the match filter
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
% preliminaries
 n=round(mlength/(t(2)-t(1)))+1;
 trin2=trin(:);
 trdsign=trdsign(:);
% generate the Toeplitz matrix for the normal equations
 TRIN= convmtx(trin2,n);
% solve the equations with left division
 if flag==1
  mfilt=TRIN\[trdsign;zeros(n-1,1)];
  tm=xcoord(0.,t(2)-t(1),mfilt);
 else
  nh=fix(n/2);
  top=[zeros(nh,1);trdsign;zeros(n-nh-1,1)];
  mfilt=TRIN\top;
  tm=xcoord(-(nh)*(t(2)-t(1)),t(2)-t(1),mfilt);
 end
 [j,k]=size(trin);
 if j==1, mfilt=mfilt.'; tm=tm'; end
   
 
 
