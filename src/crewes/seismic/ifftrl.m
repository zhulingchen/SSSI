function [r,t]= ifftrl(spec,f)
% IFFTRL inverse Fourier transform for real time series
%
% [r,t]= ifftrl(spec,f)
%
% Inverse fourier transform to a real trace. The spectrum is
% assumed to be of length floor(n/2+1) where n is the length of the time
% series to be created. (This is automatically the case if it was created by
% fftrl). This is a slightly ambiguous situation since, for
% example, time series of length 4 and length 5 will both result in 3
% frequency samples. However, only the first case will have a sample at
% Nyquist. Since the sample at Nyquist should be purely real, the algorithm
% tests the sample spec(end) for a nonzero imaginary part. If a zero
% imaginary part is found, it is assumed to be a sample at Nyquist. So, a
% problem could arise if the spectrum has been altered such that the
% Nyquist is no longer purely real. So, a good practice is to always force
% the time series to have an even number of samples (before fftrl) and then
% force the Nyquist sample to be purely real (before ifftrl).
%
% spec= input spectrum
% f= input frequency sample vector (or 2D matrix or 3D array)
% r= output trace
% t= output time vector
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

%test for matrix or vector
[m,n1,n2]=size(spec);
itr=0;
%test for vector
if( ~(m-1)*(n1-1) )
	if(m==1) spec=spec.'; itr=1; end
	nsamps=length(spec);
	nx=1;
    ny=1;
else
    %we have an array
	nsamps=m;
	nx=n1;
    ny=n2;
end

%test for 3D
%  test=size(spec);
%  if(length(test)==3)
%      threeD=1;
%  else
%      threeD=0;
%  end

% form the conjugate symmetric complex spectrum expected by ifft
%test for presence of nyquist
nyq=0;
rnyq=real(spec(end));
inyq=imag(spec(end));
small=100*eps;
if(rnyq==0 && inyq==0) nyq=1;
elseif(rnyq==0 && abs(inyq)<small) nyq=1;
elseif(rnyq==0) nyq=0;
elseif(abs(inyq/rnyq)<small) nyq=1;
end
%if(isreal(spec(end))) nyq=1; end
if(nyq)
    L1=1:nsamps;L2=nsamps-1:-1:2;
else
    L1=1:nsamps; L2=nsamps:-1:2;
end

symspec=[spec(L1,:,:);conj(spec(L2,:,:))];

% if(threeD)
%     symspec=[spec(L1,:,:);conj(spec(L2,:,:))];
% else
%     symspec=[spec(L1,:);conj(spec(L2,:))];
% end

% transform the array, force real result
 r=real(ifft(symspec));

% build the time vector
 n=size(r,1);
 df=f(2)-f(1);
 dt=1/(n*df);
 %dt=1./(2*f(length(f)));
 %tmax=size(r,1)*dt;
 %t=linspace(0.,tmax,size(r,1));
 t=dt*(0:n-1)';
 
 if(itr==1)
 	r=r.';
    t=t';
 end
