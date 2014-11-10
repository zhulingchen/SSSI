function s2p=bandwidth_xfer(s1,s2,n)
% BANDWIDTH_XFER ... transfer the bandwidth of s1 to s2
%
% s2p=bandwidth_xfer(s1,s2,n)
%
% Transfer the amplitude spectrum of s1 to s2 in a smooth way. Let A1 and
% A2 be the smoothed amplitude spectra of s1 and s2. Then the the complex
% spectrum of s2P is S2p=S2.*A1./A2, where S2 is the complex spectrum of
% S2.
%
% s1 ... first input signal
% s2 ... second input signal (must be the same length as s1)
% n ... number of samples in convolutional smoother to be applied to the
%       amplitude spectra of s1 and s2.
% ********* default is 10% of the length of the spectra ************
% s2p ... the new s2 with the amplitude spectrum of s1.
%
% 
% by G.F. Margrave, Nov. 2005
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

if(length(s1)~=length(s2))
    error('s1 and s2 must be the same length')
end

% impose the bandwidth of the first signal on the second
s1p=padpow2(s1);
s2p=padpow2(s2);
t=(0:length(s1p))*.001;%just invent a t
[S1,f]=fftrl(s1p,t);
S2=fftrl(s2p,t);
if(nargin<3)
    n=round(length(S1)/10);
end
AS1=convz(abs(S1),ones(1,n)/n);
AS2=convz(abs(S2),ones(1,n)/n);
S2p=S2.*AS1./AS2;
s2p=ifftrl(S2p,f);
s2p=s2p(1:length(s2));