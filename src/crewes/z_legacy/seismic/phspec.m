function phspec(t,s,flag,n)
% phspec(t,s,flag,n)
% phspec(t,s,flag)
% phspec(t,s)
%
% PHSPEC plots a simple Fourier phase spectrum.
%
% Note: PHSPEC automatically creates a new figure and invokes
% SIMPLEZOOM to provide interactive zooming. See help for SIMPLEZOOM
%
% s= input trace
% t= input time vector
% flag= 0... apply an n length mwindow to s before transforming
%       1... apply an n length half-mwindow to s before transforming
%       2... transform directly without windowing
% ************* default=2 **************
% n= transform length, if not supplied, then fft length
%    will be length(s)
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
figure
[m,nn]=size(s);
if(nargin<4)
   n=nn;
 end
if nargin<3
   flag=2; %window flag
end
%make sure its a row vector
if(m~=1)
	s=s.';
end
   if flag <2
      mw=ones(size(s));
      if flag==0, mw=mwindow(s(1,:));end
      if flag==1, mw=mwhalf(s(1,:));end
      for it=1:m, s(it,:)=s(it,:).*mw;end
    end
    %adjust for time zero
    izero=near(t,0);
    %make sure its close
    if(abs(t(izero))<t(2)-t(1))
    	s=[s(izero:length(s)) s(1:izero-1)];
    else
    	disp('***WARNING*** unable to find time zero, phase may be inaccurate')
    end
    
    %spectrum
    [spec,f]=fftrl(s,t);
	
	%to db & phs
	spec=todb(spec);
    %s=fft(s',n);
    
    ph=180*imag(spec)/pi;
    
    
% now plot it
   plot(f,ph);
   grid;
   set(gca,'xlabel',text(0,0,' Frequency (Hz)'));
   set(gca,'ylabel',text(0,0,' phase in degrees'),'ylim',[-180 180]);
   simplezoom 
  
