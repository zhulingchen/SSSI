function Sdb= todb(S,refamp)
% TODB: converts from (real,imaginary) to (decibels, phase)
%
% Sdb = todb(S,refamp);
% Sdb = todb(S);
%
% converts a complex spectrum from (real,imaginary) to 
% (decibels,phase angle). Thus if S is a complex spectrum,
% then plot(real(todb(S))) will plot the decibel amplitude
% and plot(imag(todb(S))) will plot the phase spectrum (radians)
%
% S= input complex spectrum
% refamp = input reference amplitude
%           decibels are dbdown from this value. If defaulted,
%           refamp=max(abs(S))
% Sdb= output complex spectrum in the special format:
%       (dbdown, phase angle)
% 
% the INVERSE of this process is the function FROMDB found in
% the seismic toolbox
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

 amp=abs(S); 
 phs=angle(S); 
if(nargin==1)
   refamp=max(max(amp));
end

amp=20*log10(amp/refamp);
Sdb=amp + i*phs;
    
