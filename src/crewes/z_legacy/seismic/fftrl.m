function [spec,f]= fftrl(s,t,percent,n)

% [spec,f]= fftrl(s,t,percent,n)
% [spec,f]= fftrl(s,t,percent)
% [spec,f]= fftrl(s,t)
%
% Forward fourier transform of a real trace. This is done in brute
% force fashion to mimic the result of Vern Herbert's real to
% complex FFT. A penalty of a factor of two is paid. Relative to
% MATLAB's fft it returns only the positive frequencies in an array
% half the size. If the input trace is a vector, then the return
% is simply the transform of that trace. If a matrix, then each 
% column of the matrix is transformed.
%
% s= input trace (or gather) 
% t= input time coordinate vector
% percent= specifies the length of a raised coisine taper to be
%          applied to s prior to any padding. Taper is a percent
%          of the length of s. Taper is applied using "mwindow"
%          from the seismic toolbox. Default=0%
% n= length to which the input trace is to be padded with zeros.
%  
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
		n=length(t);
 end
 if(nargin<3)
   percent=0.0;
 end
% determine number of traces in ensemble
 [l,m]=size(s);
 ntraces=1;
 itr=0; %transpose flag
 if(l==1) nsamps=m; itr=1; s=s(:); %switch to column vectors
 elseif(m==1) nsamps=l;
 else
 	nsamps=l; ntraces=m;
 end
 if(nsamps~=length(t))
		t=t(1)+(t(2)-t(1))*(0:nsamps-1);
		 if(nargin<4)
			n=length(t);
		 end
	%error(' time vector and trace matrix don''t match in length');
 end
 
% apply the taper
 if(percent>0)
	 mw=mwindow(nsamps,percent)';
	 mw=mw(:,ones(1,ntraces));
	 s=s.*mw;
	 clear mw;
	 %for k=1:ntraces,
		%s(:,k)=s(:,k).*mw';
	 %end
 end
% pad s if needed
 if (nsamps<n),
	s=[s;zeros(n-nsamps,ntraces)];
  	nsamps=n; 
 end

% transform the array, This used to be done in a loop to conserve memory
	%spc=ones(size(s))+i*ones(size(s));
	spec=fft(s,nsamps);
	spec=spec(1:round(n/2+1),:);% save only the positive frequencies
	clear s;
    %spec=ones(nsamps/2+1,ntraces)+i*ones(nsamps/2+1,ntraces);
    %for k=1:ntraces,
         %temp=fft(s(:,k));
         %spc(:,k)=temp;
         %spec(:,k)=temp(1:round(n/2+1)); % save only the positive frequencies
    %end
% build the frequency vector
aafnyq=find( t > 0 ); % ever heard of negative time ?
aa1=min(aafnyq); aa2=min(aafnyq)+1;
 fnyq=1./( 2*( t(aa2) - t(aa1) ) );
 nf=size(spec,1);
 f=linspace(0.,fnyq,nf)';
 
 if(itr)
 	f=f';
 	spec=spec.';
 end
