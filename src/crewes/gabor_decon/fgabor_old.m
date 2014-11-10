function [tvs,tout,fout]=fgabor_old(trin,t,twin,tinc,padflag,normflag,taperpct)
% FGABOR: forward Gabor transform with Gaussian analysis windowing
%
% [tvs,tout,fout]=fgabor_old(trin,t,twin,tinc,padflag,normflag,taperpct)
% 
% FGABOR performs a forward Gabor transform of a seismic trace using a
% Gaussian analysis window. The transform is implemented by windowing the 
% trace multiple times with a temporally shifted Gaussian window. An 
% ordinary fft is performed over each window. The output is a 2D matrix,
% called tvs, with the row coordinate being the time of the center of each
% Gaussian window and the column coordinate being the Fourier frequency in
% Hz. This tvs, or time variant spectrum, is also called the Gabor spectrum.
% The Gabor spectrum may be inverted with either IGABOR or IGABOR_SYN.
%
% trin= input trace 
% t= time coordinate vector for trin
% twin= width (seconds) of the Gaussian window
% tinc= temporal shift (seconds) between windows
% padflag= if 0, the trace is transformed without padding. If 1,
%   it is padded with zeros to the next power of 2 (unless it already is
%   a power of 2)
% ************** default = 1 ***************
% normflag= if 1, normalize the tvs for finite length time series
% ************** default = 1 ******************
% taperpct = size of taper applied to the end of the trace (max time)
%                 expressed as a percent of twin
% ************** default = 200% *********************
% tvs= output time-variant spectrum (complex valued)
% tout= column vector giving the row coordinate of tvs
% fout= row vector giving the column coordinate of tvs
%
% by G.F. Margrave, May 2001
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

if(nargin<7) taperpct=200; end
if(nargin<6) normflag=1; end
if(nargin<5) padflag=1; end

% %scale taperpct to the trace length
% taperpct=taperpct*twin/max(t);

tmin=t(1);
t=t-tmin;

%make sure we have row vectors
[m,n]=size(trin);
if(n==1) trin=trin.'; end
[m,n]=size(t);
if(n==1) t=t'; end
dt=t(2)-t(1);

% %taper trace
% if(taperpct>0)
%    trin=trin.*mwhalf(length(trin),taperpct)';
% end

%pad trace
if(padflag)
    trin=padpow2(trin);
    t=(0:length(trin)-1)*(t(2)-t(1));
end

%determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
tmax=t(end);
nwin=tmax/tinc+1; %this will generally be fractional
nwin=round(nwin);
tinc=tmax/(nwin-1); %redefine tinc

tout=(0:nwin-1)*tinc;

%loop over windows
itn=zeros(1,nwin);
if(normflag)
    normfactor=zeros(size(t));
    for k=1:nwin
        %build the gaussian
        tnot=(k-1)*tinc;
        itn(k)=round((tnot-t(1))/dt)+1;
        gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
        %if(normflag) tnorm=tnorm+gwin; end
        normfactor=normfactor+gwin;
    end
else
    normfactor=ones(size(t));
end

aaa=1;
%
gwinmat=zeros(nwin,length(t));
for k=1:nwin
    %build the gaussian
    tnot=(k-1)*tinc;
    itn(k)=round((tnot-t(1))/dt)+1;
    gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
    gwin=gwin./normfactor;
%     gwinmat(k,:)=gwin;
    %window and fft
    if(k==1)
        [tmp,fout]=fftrl(gwin.*trin,t);
        tvs=zeros(nwin,length(tmp));
        tvs(k,:)=tmp;
    elseif(k<nwin)
        tvs(k,:)=fftrl(gwin.*trin,t);
    else
        tvs(k,:)=fftrl(gwin.*trin,t);
    end
    gwinmat(k,:)=trin.*gwin;
end

tout=tout+tmin;

%normalize
% if(normflag)
% for k=1:nwin
%     tvs(k,:)=tvs(k,:)/tnorm(itn(k));
% end
% end
