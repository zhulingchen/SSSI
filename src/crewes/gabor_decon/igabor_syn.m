function [trout,t]=igabor_syn(tvs,f,twin,tinc,afilt,normflag)
% IGABOR_SYN: inverse Gabor transform with Gaussian synthesis windowing and TV filtering
%
% [trout,t]=igabor_syn(tvs,fcol,twin,tinc,afilt,normflag)
% 
% IGABOR_SYN performs an inverse Gabor transform of a Gabor spectrum. This is
% implemented with a gaussian synthesis window and an optional filter. The input 
% may be either a vector, in which case it is assumed to be an ordinary Fourier 
% spectrum, or a matrix giving
% a time-variant or Gabor spectrum. In the latter case it is immediately summed along
% the columns (collapsing the row dimension) to give an ordinary Fourier spectrum.
% The algorithm defines a series of Gaussian windows that increment regularly along
% the output time axis. With a loop over window position, the Fourier spectrum is
% multiplied by the filter corresponding to the window center, inverse transformed,
% and windowed. Then it is summed into the output vector.
%
% tvs= input time variant (Gabor) spectrum or ordinary Fourier spectrum. This is
%      typically created by FGABOR or IFFTRL.
% f= frequency coordinate vector for the columns of tvs.
% twin= width (seconds) of the Gaussian window
% tinc= temporal shift (seconds) between windows
% afilt= time-variant filter to apply. This must have the same number of columns as
%      the number of frequencies in tvs and the number of rows must equal the number
%      of gaussians.
% trout= output time series
% t= time coordinate vector for trout
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


[nwin,nf]=size(tvs);
if((nwin~=1)|(nf~=1))
    spec=sum(tvs);
else
    spec=tvs;
end
nf=length(spec);
clear tvs

%we need to do the first transform to determine the number of windows needed
[sig,t]=ifftrl(spec.*afilt(1,:),f);

%determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
tmax=t(end);
nwin=tmax/tinc+1; %this will generally be fractional
nwin=round(nwin);
tinc=tmax/(nwin-1); %redefine tinc

tnot=(0:nwin-1)*tinc; %these are the window center times

%loop over windows
%if(normflag)
%    tnorm=zeros(size(t));
%    itn=zeros(1,nwin);
%end
trout=zeros(size(sig));
for k=1:nwin
    %build the gaussian
    %tnot=(k-1)*tinc;
    %itn(k)=round((tnot(k)-t(1))/dt)+1;
    gwin=exp(-((t-tnot(k))/twin).^2)/(sqrt(pi)*twin/tinc);
    trout= trout + gwin.*sig;
    if(k<nwin)
        sig=ifftrl(spec.*afilt(k+1,:),f);
    end
end

trout=trout(:);
