function [MX,Phase,ff]=sacfft(data,tt,IFPLOT)
% function [MX,Phase,ff]=sacfft(data,tt,IFPLOT)
%
% translate the data from time domain into
% frequency domain.
% SAC Matlab package
% Zhigang Peng
% Sun Feb 25 12:49:06 PST 2001
% modified from 
% http://www.mathworks.com/support/tech-notes/1700/1702.shtml

if nargin <=2 IFPLOT = 0; end
if nargin <=1 error('Not enough input argument');
end


%if abs((t(2) - t(1)) - (t(end) - t(end-1)))>=0.000001
% error(' Freq inteval is not the same!!');
%end

delta = tt(2) - tt(1);


Fs=1/delta;                  % sampling frequency
Fn=Fs/2;                  % Nyquist frequency

% Next highest power of 2 greater than or equal to
% length(x):

NFFT=2.^(ceil(log(length(data))/log(2)));
% Take fft, padding with zeros, length(FFTX)==NFFT

FFTX=fft(data,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
% fft is symmetric, throw away second half
FFTX=FFTX(1:NumUniquePts);
MX=abs(FFTX);            % Take magnitude of X
Phase=angle(FFTX);            % Take phase of X

% Multiply by 2 to take into account the fact that we
% threw out second half of FFTX above

%MX=MX*2;

%MX(1)=MX(1)/2;
%if ~rem(NFFT,2),
%    MX(length(MX))=MX(length(MX))/2;
%end
% comment out to follow SAC and Glenn's approach

% Scale the FFT so that it is not a function of the 
% length of x.
% MX=MX/length(x);                  %matlab approach

MX=MX*delta; % Glenn's approach

ff=(0:NumUniquePts-1)'*2*Fn/NFFT;

if(IFPLOT)
  figure;
  subplot 211;
  plot(tt,data);
  xlabel('Time');
  title('Data in the Time Domain');
  subplot 212;
  loglog(ff,MX);
  xlabel('Frequency');
  title('Data in the Frequency Domain');
end
