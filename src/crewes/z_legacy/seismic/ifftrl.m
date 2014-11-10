function [r,t]= ifftrl(spec,f)
%
% [r,t]= ifftrl(spec,f)
%
% Inverse fourier transform to a real trace. This is done in brute
% force fashion to mimic the result of Vern Herbert's complex to
% real FFT. A penalty of a factor of two is paid. The spectrum is
% assumed to be 1 sample longer than a magic length (power of 2).
% This will automatically be the case if it was created by fftrl
%
% spec= input spectrum
% f= input frequency sample vector 
% r= output trace
% t= output time vector
%
% by G.F. Margrave, May 1991
%
%test for matrix or vector
[m,n]=size(spec);
itr=0;
if( ~(m-1)*(n-1) )
	if(m==1) spec=spec.'; itr=1; end
	nsamp=length(spec);
	ntr=1;
else
	nsamp=m;
	ntr=n;
end
% form the conjugate symmetric complex spectrum expected by ifft
  L1=1:nsamp;L2=nsamp-1:-1:2;
  symspec=[spec(L1,:);conj(spec(L2,:))];
% transform the array
 r=real(ifft(symspec));
% build the time vector
 dt=1./(2*f(length(f)));
 tmax=(length(r)-1)*dt;
 t=linspace(0.,tmax,length(r));
 
 if(itr==1)
 	r=r.';
 end
