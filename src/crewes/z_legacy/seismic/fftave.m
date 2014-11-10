function [spec,f]= fftave(s,t,percent,n)
% [spec,f]= fftrl(s,t,percent,n)
% [spec,f]= fftrl(s,t,percent)
% [spec,f]= fftrl(s,t)
%
% Forward fourier transform of a real trace. This is done in brute
% force fashion to mimic the result of Vern Herbert's real to
% complex FFT. A penalty of a factor of two is paid. Relative to
% MATLAB's fft it returns only the positive frequencies in an array
% half the size. FFTRL is ensemble capable.
%
% s= input trace 
% t= input time vector
% percent= specifies the length of a raised coisine taper to be
%          applied to s prior to any padding. Taper is a percent
%          of the length of s. Taper is applied using "mwindow"
%          from the seismic toolbox. Default=0%
% n= length to which the input trace is to be padded with zeros.   
% spec= output spectrum
% f= output frequency sample vector
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
 if(nargin<4)
			n=length(s);
	end
 if(nargin<3)
   percent=0.0;
	end
% determine number of traces in ensemble
 ntraces=min(size(s));
 nsamps=max(size(s));
 [l,m]=size(s);
 if l>m, s=s.'; end
% apply the taper
 if(percent>0)
  mw=mwindow(nsamps,percent);
  for k=1:ntraces,
   s(k,:)=s(k,:).*mw;
  end
 end
% pad s if needed
 if length(s)<n,
  s=[s zeros(ntraces,n-length(s))];
  nsamps=n; 
 end
% transform the array, This is done in a loop to conserve memory
 spec=zeros(1,nsamps/2+1);
 for k=1:ntraces,
     temp=fft(s(k,:));
     spec=spec + abs(temp(1:round(n/2+1)));% save only the positive frequencies
 end
% build the frequency vector
 spec = spec/ntraces;
 fnyq=1./(2*(t(2)-t(1)));
 f=linspace(0.,fnyq,length(spec));
