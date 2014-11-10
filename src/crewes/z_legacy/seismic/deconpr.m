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
% design prediction operator
 x= predict(trdsign,nop,nlag,stab);
% form the predicted part of trin
 trinhat= conv(trin,x);
% delay and subtract
 trout= trin(nlag+1:length(trin))-trinhat(1:length(trin)-nlag);
 trout= [trin(1:nlag) trout];
% balance the output
 trout=balans(trout,trin);
