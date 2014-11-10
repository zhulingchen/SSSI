function [r,r2,w,q,q2]= deconz(trace,trace_design,nw,stab)
%
%
% [r,w,q]=deconz(trace,trace_design,nw,nr,stab)
%
%
% routine performs a Ziolkowski style deconvolution of the input trace
%
% trace= input trace to be deconvolved
% trace_design= input trace to be used for operator design
% nw= expected temporal length of the wavelet estimate. This is equivalent
%       specifying the number of autocorrelation lags in deconw
% ************** default = round(length(trace_design)/10) **************
% stab= stabilization factor expressed as a fraction of the
%       zero lag of the autocorrelation.
%      ********* default= .0001 **********
%
% r= estimate reflectivity (solution by left division of conv mtx)
% r2= estimate reflectivity (solution by normal equations)
% w= wavelet as estimated by Wiener decon
% q= Ziolkowski's quality factor for r
% q2= Ziolkowski's quality factor for r2
%
% by: G.F. Margrave, May 1991
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
  if nargin<4
    stab=.0001;
  end
  if(nargin<3)
      nw=round(length(trace_design)/10);
  end
%step 1 is to estimate the Wavelet in the Wiener style
% generate the autocorrelation
  a=auto(trace_design,nw,0);
% hanning window applied to auto
  window=hanning(2*length(a)-1);
  a=a.*window(length(a):end)';
% stabilize the auto
  a(1)=a(1)*(1.0 +stab);
  %a=a/a(1);%normalize auto
% generate the right hand side of the normal equations
  b=[1.0 zeros(1,length(a)-1)];
% do the levinson recursion
  x=levrec(a,b);%so x is the inverse of the desired wavelet
% normalize the inverse operator
  %x=x/sqrt(x'*x);
% now, invert x to get the wavelet estimate
  w=ifft(1./fft(x));
%
% Step 2: the reflectivity is estimated by a match filter process where 
% sum((convm(w,r)-trace).^2) is minimized
% 

% the length of r is predetermined by the convolution rule for w*r=s
% that is length(s)=length(w)+length(r)-1
%
  nr=length(trace)-length(w)+1;
  
%make a convolution matrix from w
  W=convmtx(w,nr);
  r=pinv(W,.1)*trace;

%form the normal equations
A=W'*W;%autocorrelation matrix for wavelet
D=W'*trace;%crosscorrlelate trace with wavelet
%stabilize A
A=A+diag(diag(A)*(1+stab/10));

%solve
r2=A\D;
 
  % quality factor
  trace_ideal=conv(w,r);
  trace_ideal2=conv(w,r2);
  q=sum(trace_ideal.^2)/sum(trace.^2);
  q2=sum(trace_ideal2.^2)/sum(trace.^2);
  

  
