function wnorm = wavenorm(w,tw,flag)
% WAVENORM: normalize a wavelet
%
% wnorm = wavenorm(w,tw,flag)
%
% WAVENORM normalizes a wavelet by one of several criteria.
% The choices are: (1) normalize the maximum absolut value to unity
% (2) normalize such that a sine wave at the dominant frequency passes
% with unit amplitude (3) normalize the rms amplitude to unity.
%
% w ... input wavelet
% tw ... time coordinate vector for w
% flag ... (1) normalize the maximum absolute value to unity
%          (2) normalize such that a sine wave at the dominant frequency 
%              passes with unit amplitude
%          (3) normalize the rms amplitude to unity
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

if(flag==1)
    wnorm = w/max(abs(w));
elseif(flag==2)
    %compute spectrum
    [W,f]=fftrl(w,tw);
    A=real(todb(W));
    %pick maximum amplitude
    ind=find(A==max(A));
    fdom = f(ind(1));
    refwave=sin(2*pi*fdom*tw);
    reftest=convz(refwave,w);
    fact=max(refwave)/max(reftest);
    wnorm=w*fact;
elseif(flag==3)
    rms_ave=norm(w)/sqrt(length(w));
    wnorm=w/rms_ave;
else
    error('invalid flag')
end
