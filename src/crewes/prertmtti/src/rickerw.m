function [w,tw] = rickerw(fpeak,dt)
% rickerw(fpeak,dt)
%
% rickerw: ricker wavelet
%
% 	fpeak:    the peak frequency
%
%	dt:       time sampling
%
%	w:        ricker wavelet
%
%	tw:       time sequence
%
%   Xiang Du, May   2007, $1.0
%             Sept. 2008, $1.1
%             Nov.  2008, $1.2

 nw=2.2/fpeak/dt;
 nw=2*floor(nw/2)+1;
 nc=floor(nw/2);
 w = zeros(nw,1);

 k=[1:1:nw]';

 alpha = (nc-k+1).*fpeak*dt*pi;
 beta=alpha.^2;
 w = (1.-beta.*2).*exp(-beta);

if nargout>1;
  tw = -(nc+1-[1:1:nw])*dt;
end
