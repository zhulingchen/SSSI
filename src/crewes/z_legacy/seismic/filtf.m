 function trout=filtf(trin,t,fmin,fmax,phase,max_atten)
%
% trout=filtf(trin,t,fmin,fmax,phase,max_atten)
% trout=filtf(trin,t,fmin,fmax,phase)
% trout=filtf(trin,t,fmin,fmax)
%
% FILTF filters the input trace in the frequency domain.
% Trin is automatically padded to the next larger power of
% two and the pad is removed when passing trout to output. 
% Filter slopes are formed from Gaussian functions. (If the
% input is already a power of 2, no pad is applied. No other 
% protection against filter wrap around is provided)
%
% trin= input trace
% t= input trace time coordinate vector
% fmin = a two element vector specifying:
%        fmin(1) : 3db down point of filter on low end (Hz)
%        fmin(2) : gaussian width on low end
%   note: if only one element is given, then fmin(2) defaults
%         to 5 Hz. Set to [0 0] for a low pass filter  
% fmax = a two element vector specifying:
%        fmax(1) : 3db down point of filter on high end (Hz)
%        fmax(2) : gaussian width on high end
%   note: if only one element is given, then fmax(2) defaults
%         to 10% of Fnyquist. Set to [0 0] for a high pass filter. 
% phase= 0 ... zero phase filter
%       1 ... minimum phase filter
%  ****** default = 0 ********
% note: because the Fourier transform pad is removed by simply 
%	truncating the filtered trace, there is an apparent limit 
%	placed on the attenuation of the rejected spectrum. This is
%	easily circumvented by applying an appropriate window after
%	filtering. To see this effect compare:
% 	dbspec(t,trout) and dbspec(t,trout.*mwindow(trout))
% note: Minimum phase filters are approximate in the sense that
%  the output from FILTF is truncated to be the same length as the
%  input. This works fine as long as the trace being filtered is
%  long compared to the impulse response of your filter. Be wary
%  of narrow band minimum phase filters on short time series. The
%  result may not be minimum phase.
% 
% max_atten= maximum attenuation in decibels
%   ******* default= 80db *********
%
% trout= output trace
%
% by G.F. Margrave, May 1991
% 
 
% set defaults
 if nargin < 6
   max_atten=80.;
 end
 if nargin < 5
   phase=0;
 end
 if length(fmax)==1
   fmax(2)=.1/(2.*(t(2)-t(1)));
 end
 if length(fmin)==1
   fmin(2)=5;
 end
% convert to column vectors
[rr,cc]=size(trin);
trflag=0;
if(cc>1)
		trin=trin';
		t=t';
		trflag=1;
end
% forward transform the trace
  ntrout=length(trin);
  trin=padpow2(trin);
  t=xcoord(t(1),t(2)-t(1),length(trin));
  [Trin,f]=fftrl(trin,t);
% design filter
[fltr,f]=filtspec(t(2)-t(1),t(length(t))-t(1),fmin,...
		fmax,phase,max_atten);
% apply filter
  trout=ifftrl( Trin.*fltr(1:length(f)),f );
  trout=trout(1:ntrout);
  if(trflag)
		trout=trout';
  end
