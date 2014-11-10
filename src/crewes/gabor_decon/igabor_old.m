function [trout,t]=igabor_old(tvs,fcol)
% IGABOR: inverse Gabor transform without synthesis windowing
%
% [trout,t]=igabor_old(tvs,fcol)
% 
% IGABOR performs an inverse Gabor transform of a Gabor spectrum. This is
% implemented with a synthesis window of unity. For a gaussian
% synthesis window use IGABOR_SYN. The algorithm is a simple sum
% along each column (collapsing the row dimension) of the Gabor
% spectrum followed by an ordinary IFFT.
%
% tvs= input time variant spectrum or Gabor spectrum. This is typically
%      created by FGABOR.
% fcol= frequency coordinate vector for the columns of tvs
% trout= output time series
% t= time coordinate vector for trout
%
% by G.F. Margrave, May 2001
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


[nwin,nf]=size(tvs);

%loop over windows
%for k=1:nwin
%    if(k==1)
%        [trout,t]=ifftrl(tvs(k,:),fcol);
%    else
%        trout=trout+ifftrl(tvs(k,:),fcol);
%    end
%end

% the one liner achieves the same thing as the loop because the fft and sum
% commute
[trout,t]=ifftrl(sum(tvs),fcol);

%trout=trout(:)/nwin;
trout=trout(:);
