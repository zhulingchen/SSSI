function eqndefinit(transfer,hmasterfig,lognames,logdefaults,coefs,...
		exps,newlogdes)
% eqndefinit(transfer,hmasterfig,lognames,logdefaults,coefs,...
%		exps,newlogdes)
%
% EQNDEFINIT is called by LOGEDIT to initiate a dialog  to define the
% computation of a new log from existing ones via an algebraic equation.
% The equation is of the form ((a*l1^k + b*l2^m)/(c*l3^n+d*l4^o))^p 
% Here l1,l2,l3,l4 are logs and the other terms are scalar constants. The
% dialog allows the definition of the logs, the scalars, and the destination of
% of the computed log.
%
% transfer ... transfer command to be called when the user terminates
%		the dialog
% hmasterfig ... handle of the masterfigure (usually a LOGEDIT window) in
%		control of the dialog
% lognames ... list of names of logs with individual names separated by '|' 
%		such as: 'fred|sam|billy|wilma'
% logdefaults ... vector of length four containing four integers which are the 
% 		defaults for the four logs. These are taken to be row indicies into
%		lognames
% coefs ... vector of length four containing the defaults for a,b,c, and d
% exps ... vector of length five containing the defaults for k,m,n,o, and p
% newlogdes ... default descriptor for the new log. (Its "name" will
%	be automatically implied by the mnemonic form its type.)
%
%
% G.F. Margrave, Jan 1996
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
% pack the information into the current axes userdata
hax=get(hmasterfig,'currentaxes');
figure(hmasterfig);
set(hax,'userdata',[abs(transfer) nan abs(lognames) nan hmasterfig,logdefaults,...
	coefs,exps nan abs(newlogdes)]);
eqndef('init')
