function [mask,htrin]=hmask(trin,thresh)

% [mask,htrin]=hmask(trin,thresh)
% [mask,htrin]=hmask(trin)
%
% Compute the Hilbert Reflectivity Mask of trin.
% HMASK returns a trace which is zero everywhere except at the
% samples corresponding to peaks of abs(hilbert(trin)) where it
% is 1.0
%
% trin= input trace
% mask= output hilbert mask trace
% htrin= complex hilbert transform of trin. (e.g. hilbert(trin))
% thresh = threshold significance for peaks. Peaks on the Hilbert
%	envelope are ignored if they are smaller than thresh times
%	the maximum peak.
%	********** default = .05 *********
% 
% by G.F. Margrave, July 1991
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

% $Id: hmask.m,v 1.4 2004/07/30 22:03:45 kwhall Exp $

if(nargin<2)
		thresh=.05;
end

% compute hilbert transform and envelope
htrin=hilbert(trin).';
env=abs(htrin);
% build hilbert mask
iex=findex(env);
 mask=zeros(size(env));
 if length(iex)>=1,
  mask(iex)=ones(size(iex));
 end
 mm=max(env);
 iz=find(env<mm*thresh);
 mask(iz)=zeros(size(iz));
