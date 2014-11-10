function [lognums,coefs,exps,type,des]=eqndeffini

% [lognums,coefs,exps,type,des]=eqndeffini
%
% Call this to complete the Log algebra equation definition dialog. 
% Return values mean:
%  lognums = vector of length four giving the numbers of the four logs
%  coefs = vector of length four giving the scalar coeficients of the four logs
%  exps = vector of length five giving the scalar exponents of the four 
%         logs and the overall exponent
% type = type flag for the new log
% des = string descriptor for the new log
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

	params=get(gca,'userdata');
	if( params==-1 ) %test for a cancel
		lognums=-1;coefs=[];exps=[];type=[];des=[];
		return;
	end

	lognums=params(1:4);
	coefs=params(5:8);
	exps=params(9:13);
	type=params(14);
	
	des= setstr(params(15:length(params)));
