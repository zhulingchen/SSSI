function freqloc=gaborfreq(s,t,twin,tinc,p)
% GABORFREQ ... calculates a local dominant frequency using the Gabor transform
%
% freqloc=gaborfreq(s,t,twin,tinc,p)
%
% An input signal, s, is Gabor transformed using fgabor.  This gives a
% time-frequency decomposition which has a local Fourier spectrum at each
% of the times t0=0, tinc, 2*tinc, ... max(t). The spectral resolution of
% these is determined by twin which is the half-width of a temporal
% Gaussian window. The dominant frequency of each of these local spectra is
% then calculated by fdom=integral(f.*spec.^p)/integral(spec.^p). In this
% expression, "integral" is an integration over all frequencies, f is
% frequency, and p is a small integrer usually 1 or 2. (fdom is often
% called the "centroid frequency".) The result is a local frequency which
% is assigned to each window center. As a final step, these local
% frequecies are interpolated to the time coordinate of the input signal s.
%
% s... input signal
% t... time coordinate vector for s in seconds
% twin ... Gaussian window half width in seconds. A good value here might
%       be 2-10% of the signal length.
% tinc ... increment between Gaussians in seconds. Typically this is
%       smaller than twin by a factor of 2-10.
% p ... spectral exponent to use in calculating the local frequency in each
%       Gabor window.
%   ****** default 2 ****** (Means use Power spectrum)
%
% freqloc ... local frequency computed from Gabor slices. freqloc will be
%       the same size as s and can be plotted by plot(t,freqloc)
%
% by G.F. Margrave, 2013
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

[S,trow,fcol]=fgabor(s,t,twin,tinc);

fdom=zeros(size(trow));
for k=1:length(trow)
      fdom(k)=sum(fcol.*abs(S(k,:)).^p)/sum(abs(S(k,:)).^p);
end

%interpolate
freqloc=interpextrap(trow,fdom,t,0);