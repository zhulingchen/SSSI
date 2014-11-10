function [slices,trow,tcol,gg]=gaborslice(signal,t,twin,tinc,p,gdb,normflag,pow2option)
% GABORSLICE: produces a set of Gabor slices using modified Gaussian windows
%
% [slices,trow,tcol,gg]=gaborslice(signal,t,twin,tinc,p,gdb,normflag,pow2option)
% 
% A Gabor slice is the result of multiplying a signal by a Gaussian whose
% center is at some particular time. A complete set of slices can be
% constructed such that they sum to recreated the signal. In order for the
% recreatiion to be exact, the Gaussians must be slightly modified (see
% gaussian_upou). The set of modified Gaussians sum exactly to one over the
% length of the input signal and are said to form a partition of unity
% (POU). An fft over the rows of the output gives a forward Gabor
% transform.
%
% signal= input trace 
% t= time coordinate vector for signal
% twin= width (seconds) of the Gaussian window
% tinc= temporal increment between windows
% p= exponent to define analysis windows. If g is a modified Gaussian
%   forming a POU, then the analysis window is g^p
%  ************ default p=1 **************
% gdb= decibel cutoff for Gaussian slices. See gaussian_upou for a
%   description. Set to inf for slices the same length as the input
%   signal.
% ************ default = inf *************
% normflag= set to zero for an approximate POU without normaization. Set to
%   1 for an exact POU. The normalization required to produce and exact POU
%   can distort amplitudes at the ends of the signal and this may not be
%   always acceptable.
% ************ default =0 ****************
% pow2option= 1 means expand signal length to a power of two by zero
%   padding, 0 means do not do this.
% ************ default = 0 ***************
% slices= matrix of Gabor slice slices, one slice per column. The number of
%   rows is the same as the number of samples in each slice, while the number
%   of columns is determined by the number of windows in the POU.
% trow= time coordinate for slice (a single column of slices)
% tcol= time coordinate of each slice (the gaussian window center time)
% gg= matrix the same size as slices containing the gaussian windows in
%   each column used to create the slices.
%
% by G.F. Margrave, July 2009
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

if(nargin<5); p=1; end
if(nargin<6); gdb=inf; end
if(nargin<7); normflag=0; end
if(nargin<8); pow2option=0; end
if(gdb<0)
    gdb=abs(gdb);
end
tmin=t(1);
%make sure we have row vectors
[m,n]=size(signal);
if(n==1); signal=signal.'; end
[m,n]=size(t);
if(n==1); t=t'; end
%build first window and POU norm factor
if(normflag==0)
    norm_factor=ones(size(t));
else
    norm_factor=0;
end
[gwin,norm_factor,tnotvec,nwin,i1]=gaussian_upou(t,tmin,twin,tinc,norm_factor,0,gdb,pow2option);

if(p==0)
    gwin=ones(size(gwin));
elseif(p<1)
    gwin=gwin.^p;
end
% 
% tcol=((0:nwin-1)*tinc+tmin)';
tcol=tnotvec;

trow=t;
%loop over windows
nt=max([length(signal) length(gwin)]);
slices=zeros(nt,nwin);
gg=slices;
for k=1:nwin
    %window
    gg(:,k)=gwin;
    slices(i1:i1+length(gwin)-1,k)=(gwin.*signal);
    if(k<nwin)
        if(p==0)
            gwin=ones(size(gwin));
        else
            %build the next gaussian
            gwin=gaussian_upou(t,tnotvec(k+1),twin,tinc,norm_factor,tnotvec,gdb,pow2option);
            if(p~=1)
                gwin=gwin.^p;
            end
        end
    end
end
