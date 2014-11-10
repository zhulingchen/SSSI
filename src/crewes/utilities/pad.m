function trout=pad(trin,trdsign,flag)
% PAD: pads (truncates) one trace with zeros to be the length of another
%
% trout=pad(trin,trdsign,flag)
%
% PAD pads (or truncates) trin to the same length as design vector trdsign
%
% trin= input trace to be padded (truncated)
% trdsign= design trace to give desired length (always a vector)
% flag=0 the pad is added to the end of trin
%     =1 the pad is added such that the central sample of trin
%        stays in the middle.
% ***********default=0 *************
% trout= output trace
%
% Example
%  pad([1 2 3 4 5], zeros(11), 1)
%  ans =
%     0     0     0     1     2     3     4     5     0     0     0
%
%
% by G.F. Margrave, June 1991
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

 if nargin<3, flag=0; end
 [m,n]=size(trin);
 if((m-1)*(n-1)) error('pad only works for vectors, not matrices'); end
 trflag=0;
 if(m==1) trin=trin.'; trflag=1; end
 nin=length(trin);
 nout=length(trdsign);
 
 %
 % we have different cases for odd and even for both input and
 % output sizes
 %     
 trout=zeros(nout,1);
 if(flag==0)
     %simple case, pad on end
     if(nout>=nin)
         trout(1:nin)=trin;
     else
         trout=trin(1:nout);
     end
 else
     %harder case. Maintain position of center sample
     %
     % we have different cases for odd and even for both input and
     % output sizes
     %
     if(iseven(nin))
         if(iseven(nout))
             %even in even out
             if(nout>nin)
                 %sample at nin/2+1 goes to nout/2+1
                 j0=(nout-nin)/2;
                 trout(j0+1:j0+nin)=trin;
             else
                 %as before
                 k0=(nin-nout)/2;
                 trout=trin(k0+1:k0+nout);
             end
         else
             %even in odd out
             if(nout>nin)
                 %sample at nin/2 goes to (nout+1)/2
                 j0=(nout+1)/2-nin/2-1;
                 trout(j0+1:j0+nin)=trin;
             else
                k0=-(nout+1)/2+nin/2+1;
                trout=trin(k0+1:k0+nout);
             end
             
         end
     else
         if(iseven(nout))
             %odd in even out
             if(nout>nin)
                 %sample at (nin+1)/2 goes to nout/2+1
                 j0=nout/2+1-(nin+1)/2;
                 trout(j0+1:j0+nin)=trin;
             else
                 %as before
                 k0=-nout/2-1+(nin+1)/2;
                 trout=trin(k0+1:k0+nout);
             end
         else
             %odd in odd out
             if(nout>nin)
                 %sample at nin/2 goes to (nout+1)/2
                 j0=(nout-nin)/2;
                 trout(j0+1:j0+nin)=trin;
             else
                k0=(nin-nout)/2;
                trout=trin(k0+1:k0+nout);
             end
             
         end
         
     end
 end
 
 if(trflag)
     trout=trout.';
 end
 
 function flag=iseven(n)
 %test for even/odd
 if(floor(n/2)*2==n)
     flag=1;
 else
     flag=0;
 end
         
     








