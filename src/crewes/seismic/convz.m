function s=convz(r,w,nzero,nout,flag)
% CONVZ: convolution then truncation for non-min phase filters
%
% s= convz(r,w,nzero,nout,flag)
% s= convz(r,w,nzero,nout)
% s= convz(r,w,nzero)
% s= convz(r,w)
%
% CONVZ is designed for a convenient convolution of a seismic
% trace with a zero phase (no time delay) wavelet. This can 
% actually be used with any non-causal wavelet by specifying the
% nzero sample which is the sample number of zero time. It
% defaults to the middle of w but can be placed anywhere. 
% If the first input argument is either a row or column vector, then 
% the output will be a similar vector. If the first argument is a matrix,
% then the output is a matrix of similar size where w has been convolved
% with each column of r.
% Uses MATLAB's CONV function.
%
% s= output trace of length nout
% r= input trace (reflectivity)
% w= input wavelet
% nzero= sample number of zero time value for wavelet
%  *********** default=round((length(wavelet)+1)/2)
% nout= length of output trace. 
%   ********** default=length(r)
% flag= 1 --> apply a cosine taper at the beginning and end of the
%            output trace
%     = 0 --> don't apply a taper
%      ********* default= 0 **********
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


 [nsamps,ntr]=size(r);
 %convert to column vectors
 transpose=0;
 if(nsamps==1); r=r.'; ntr=1; nsamps=length(r); transpose=1;end
 if(nargin<=4)
     flag=0;
 end
 if(nargin<=3)
     nout=nsamps;
 end
 if(nargin==2)
     nzero=round((length(w)+1)/2);
 end
%
% 
w=w(:);
s=zeros(nout,ntr);
for k=1:ntr
    temp=conv(r(:,k),w);
    if(flag~=1)
        s(:,k)=temp(nzero:nout+nzero-1);
    else
        s(:,k)=temp(nzero:nout+nzero-1).*(mwindow(nout,4).');
    end
end

if(transpose)
	s=s.';
end
   
