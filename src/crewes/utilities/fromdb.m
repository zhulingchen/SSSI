function S= fromdb(Sdb,refamp)
% FROMDB: convert from (db,phase) to (real,imaginary)
%
% S = fromdb(Sdb,refamp);
% S = fromdb(Sdb);
%
% FROMDB converts a complex spectrum from (decibels,phase angle) to 
% (real,imaginary). Thus if Sdb is a complex spectrum created by
% TODB, then the original complex spectrum can be recreated
% by fromdb(Sdb,refamp) if refamp is known or, to within a scale
% factor if refamp is not known.
%
% Sdb= input complex spectrum in (dbdown, phase angle) format
% refamp = input reference amplitude
%           decibels are dbdown from this value. If defaulted,
%           refamp= 1.0
% S= output complex spectrum in the (real, imaginary) format
%     
% 
% the INVERSE of this process is the function TODB found in
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

if(nargin==1)
   refamp=1.0;
 end
ten=10.0;
amp= refamp*ten.^(real(Sdb)/20.);
S= amp.*cos(imag(Sdb)) + i*amp.*sin(imag(Sdb));
