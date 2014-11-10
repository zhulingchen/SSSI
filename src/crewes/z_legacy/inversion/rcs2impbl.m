function imp=rcs2impbl(rcs,t,fmin,fmax,delflow,delfhigh)
% imp=rcs2impbl(rcs,t,fmin,fmax,delflow,delfhigh)
%
% RCS2IMPBL takes a reflectivity vector specified in time and computes
% a corresponding impedence by integrating and exponentiating it.
% The integration is bandlimited (done by BLINT) so as not to generate
% significant power outside the band [fmin fmax]. The resultant impedance
% estimate is returned with mean removed.
%
% rcs... column vector of reflection coeficients
% t  ... column vector of times (note that rcs must be regularly sampled)
% znot ... initial impedance (at t=t(1)) estimate
% fmin ... minimum frequency to pass
% fmax ... maximum frequency to pass
% delflow ... gaussian rolloff width on low end
%         ******* default min([5 flow/5]) *******
% delfhigh ... gaussian rolloff width on high end
%         ******* default min([(fnyquist-fhigh)/10 10]) *******
%
% G.F. Margrave May 1995
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
 %integrate the rcs
 if(nargin==6)
 	ircs=blint(rcs,t,fmin,fmax,delflow,delfhigh);
 elseif(nargin==5)
 	ircs=blint(rcs,t,fmin,fmax,delflow);
 else
 	ircs=blint(rcs,t,fmin,fmax);
 end
 %exponentiate
 imp= exp(2*ircs );
 imp=imp-mean(imp);
