function [trout,pefilt]= deconb(trin,trdsign,l)
%
% [trout,pefilt]=deconb(trin,trdsign,l)
%
% routine performs a Burg scheme deconvolution of the
% input trace
%
% trin= input trace to be deconvolved
% trdsign= input trace to be used for operator design
% l= prediction error filter length (and length of
%    inverse operator
%
% trout= output trace which is the deconvolution of trin
% pefilt= output inverse operator used to deconvolve trin
%
% by: G.F. Margrave, May 1991
trflag=0;
[irow,icol]=size(trin);
if(icol==1)
		trin=trin';
		trdsign=trdsign';
		trflag=1;
	end
  
% generate the prediction error filter
  pefilt=burgpr(trdsign,l);
% convolve the prediction error filter with the input trace
  trout=convm(trin,pefilt);
  trout=balans(trout,trin);
 if(trflag)
	trout=trout';
	end
