function coef=cos_taper(sp,ep,samp)
% COS_TAPER: used by KIRK_MIG
%  coefficients of cos taper
%
%  coef=cos_taper(sp,ep,samp)
%
%  coef:  coefficient with length (start-end)/samp + 1
%  sp:    start point to begin taper
%  ep:    end point to stop taper
%  samp:  the sample interval; default is 1.
%
%  By Xinxiang Li, CREWES project, U of C.
%  FEB. 1996
%
dd=[];
if nargin < 2 error('At least two input arguments needed!'); end
if nargin < 3 samp = 1. ; end
if samp < 0 samp = -samp ; end
len = abs(ep-sp)/samp;
len = len+1;
if len <= 1   coef = [1.0]; end
if len > 1
   coef=(1:len)*0.;
   dd = 1.0/(len-1)*pi*0.5;
   for i = 1:len
       coef(i) = cos((i-1)*dd);
   end
end
clear len,dd;
    
