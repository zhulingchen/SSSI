function trout=filtt(trin,t,fmin,fmax,phase,npts)
% trout=filtt(trin,t,fmin,fmax,phase,npts)
% trout=filtt(trin,t,fmin,fmax,phase)
% trout=filtt(trin,t,fmin,fmax)
%
% FILTT filters the input trace in the time domainm using a
% finite impulse response filter designed using FIR1 from the
% signal toolbox. If a minimum phase filter is called for
% it is computed from its zero phase equivalent using TOMIN from
% the seismic toolbox.
%
% trin= input trace
% t= input trace time coordinate vector
% fmin = the lowend filter cutoff
% fmax = the highend filter cutoff
% phase= 0 ... zero phase filter
%       1 ... minimum phase filter
%  ****** default = 0 ********
% npts= number of points in the filter
%   ******* default= 64 *********
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
 
% set defaults
 if nargin < 6
   npts=64;
 end
 if nargin < 5
   phase=0;
 end
 if length(fmax)>1
   fmax=fmax(1);
 end
 if length(fmin)>1
   fmin=fmin(1);
 end
% make filter
  fnyq=1./(2.*(t(2)-t(1)));
  fltr=fir1(npts-1,[fmin/fnyq fmax/fnyq]); % problem here 
  if phase==1
    fltr=tomin(fltr);
  end
% apply filter
  if phase==0
    trout=convz(trin,fltr);
    trout=trout.*mwindow(length(trout),10);
  else
    trout=convm(trin,fltr);
  end  
  
