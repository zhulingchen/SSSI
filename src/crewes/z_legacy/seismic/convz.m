function s=convz(r,w,nzero,nout,flag)

% s= convz(r,w,nzero,nout,flag)
% s= convz(r,w,nzero,nout)
% s= convz(r,w,nzero)
% s= convz(r,w)
%
% function is designed for a convenient convolution of a seismic
% trace with a zero phase (no time delay) wavelet. This can 
% actually be used with any non-causal wavelet by specifying the
% nzero sample which is the sample number of zero time. It
% defaults to the middle of w but can be placed anywhere. Also, 
% this is designed to produce an output vector of length equal
% to the first input vector (r). Uses MATLAB's CONV function.
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

 if(nargin<=4)
     flag=0;
  end
 if(nargin<=3)
     nout=length(r);
 end
 if(nargin==2)
     nzero=round((length(w))/2);
 end
%
% 
%convert to column vectors
[a,b]=size(r);
if(a==1) r=r.'; end
w=w(:);

temp=conv(r,w);
s=temp(nzero:nout+nzero-1);
if(flag==1)
   s=s.*(mwindow(nout,4).');
end

if(a==1)
	s=s.';
end
   
