function f=freqfft(t,n,flag) 
 
% f=freqfft(t,n,flag) 
% f=freqfft(t,[],flag)
% f=freqfft(t);
%
% returns the appropriate frequency coordinate vector for a 
% time series, s, whose spectrum, S, was calculated with
% S= fftshift(fft(s,n)) (sort of, see flag parameter)
%
% t= input time coordinate vector corresponding to s
% n= padded length of fft 
% ******** default n=length(t) ***********
% flag ... if 0, the f vector is returned unwrapped (f=0 in center)
%          if 1, the f vector is returned wrapped (f=0 is first sample)
%   ******* default 0 *******
% f= output frequency vector 
%
% by G.F. Margrave, 2003-2005
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

% set default
 if nargin==1
  n=length(t);
 end
 if(isempty(n))
     n=length(t);
 end
%
if(nargin<3)
    flag=0;
end
 %fnyq=.5/(t(2)-t(1));
 %df=2.*fnyq/n;
 dt=t(2)-t(1);
 df=1/(n*dt);
 if(floor(n/2)*2==n)
    %even case
    f1=(0:n/2)*df;%positive f's
    f=[-fliplr(f1(2:end)) f1(1:end-1)];
else
    %odd case
    f1=df*(0:(n-1)/2); %positive f's
    f=[-fliplr(f1(2:end)) f1];
end
if(flag)
    ind=find(f==0);
    f=[f(ind:end) f(1:ind-1)];
end