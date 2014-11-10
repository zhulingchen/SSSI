function sect = sectfilt(sect,t,fmin,fmax,phase,max_atten)
% SECTFILT: applies a bandpass filter to a seismic section
%
% sectout = sectfilt(sectin,t,fmin,fmax,phase,max_atten)
%
% SECTFILT applies a bandpass filter to a seismic section.
% It simply loops over columns and calls filtf for each trace.
%
% sectin ... input section of size nsamp x ntr. That is one trace per
%	column.
% t ... nsamp long time coordinate vector for sectin
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
% sectout ... output section of size nsamp x ntr.
%
% G.F. Margrave, CREWES Project, University of Calgary, 1996

[nsamp,ntr]=size(sect);

if(nargin < 5)
	for k=1:ntr
		sect(:,k) = filtf( sect(:,k),t,fmin,fmax);
	end
elseif( nargin < 6)
	for k=1:ntr
		sect(:,k) = filtf( sect(:,k),t,fmin,fmax,phase);
	end
else
	for k=1:ntr
		sect(:,k) = filtf( sect(:,k),t,fmin,fmax,phase,max_atten);
	end
end
