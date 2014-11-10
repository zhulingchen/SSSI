function [trout,x]= deconpr(trin,trdsign,nop,nlag,stab)
% [trout,x]= deconpr(trin,trdsign,nop,nlag,stab)
% [trout,x]= deconpr(trin,trdsign,nop,nlag)
%
% DECONPR performs Wiener predictive deconvolution by calling
% PREDICT to design a prediction filter, nop long with lag nlag
% and stab factor, using trdsign. The predicted part of trin, trinhat,   
% then formed by convolving the prediction operator with trin,
% and trout is computed by delaying trinhat by nlag samples and 
% subtracting it from trin. The prediction operator is returned
% in x.
%
% trin= input trace to be deconvolved
% trdsign= input trace used to design the prediction operator
% nop= number of points in the prediction operator
% nlag= prediction lag distance in samples
% stab= stabilazation factor expressed as a fraction of the zero
%       lag of the autocorrelation.
%  ************ default= .0001 ***********
%
% trout= deconvolved output trace
% x= prediction operator
%
% NOTE: if trin is a matrix, then each column is deconvolved and trdsign
%   must have the same number of columns as trin.
%
% See also: Peacock and Treitel, Geophysics vol 34, 1968
%  and the description of PREDICT
%
% by G.F. Margrave, July 1991
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
 if nargin<5
   stab=.0001;
 end
%test for row vector input
[m,n]=size(trin);
if((m-1)*(n-1)==0)
    %single trace
    % design prediction operator
     x= predict(trdsign,nop,nlag,stab);
     %x= predict_alt(trdsign,nop,nlag,stab);
    % form the predicted part of trin
     trinhat= conv(trin,x);
    % delay and subtract
     trout= trin(nlag+1:length(trin))-trinhat(1:length(trin)-nlag);
     if(n==1)
         trout= [trin(1:nlag);trout];
     else
         trout= [trin(1:nlag) trout];
     end
     %trout=balans(trout,trin);
else
    %a gather, we deconvolve the columns
    [mm,nn]=size(trdsign);
    if(n~=nn)
        error(' for multichannel input, the design array must have the same number of columns as the trace array')
    end
    if(length(nop)~=n)
        if(length(nop)~=1)
            error('nop must either be one value or one value per trace')
        end
        nop=nop*ones(1,n);
    end
    if(length(nlag)~=n)
        if(length(nlag)~=1)
            error('nlag must either be one value or one value per trace')
        end
        nlag=nlag*ones(1,n);
    end
    if(length(stab)~=n)
        if(length(stab)~=1)
            error('stab must either be one value or one value per trace')
        end
        stab=stab*ones(1,n);
    end
    trout=zeros(size(trin));
    for k=1:n
        x= predict(trdsign(:,k),nop(k),nlag(k),stab(k));
        %x= predict_alt(trdsign(:,k),nop(k),nlag(k),stab(k));
        % form the predicted part of trin
         trinhat= conv(trin(:,k),x);
        % delay and subtract
         tmp= trin(nlag(k)+1:length(trin),k)-trinhat(1:length(trin)-nlag(k));
         trout(:,k)= [trin(1:nlag(k),k);tmp];
         
    end
end
