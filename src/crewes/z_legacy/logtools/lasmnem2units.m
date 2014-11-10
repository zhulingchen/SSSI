function units=lasmnem2units(mnem,metric)
% units=lasmnem2units(mnem,metric)
%
% Given a LAS log mnemonic return a string giving the standard
% units for that log to be in. Set metric to 1 for metric units,
% 0 for imperial units, and -1 for time. Algorithm currently
% returns 'UNKN' for all time log units
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
if(length(mnem)<4) mnem=[mnem blanks(4-length(mnem))]; end
mnem=upper(mnem);
%strip off any number from the end of the mnemonic
id=str2num(mnem(4));
if(~isempty(id))
		mnem=mnem(1:3);
end
lm=length(mnem);
units='UNKN';
if(metric==1)
	test='DEPT';
	if(strcmp(mnem,test(1:lm))) units='M   '; return; end
	test='AU  ';
	if(strcmp(mnem,test(1:lm))) units='US/M'; return; end
	test='DT  ';
	if(strcmp(mnem,test(1:lm))) units='US/M'; return; end
	test='SON ';
	if(strcmp(mnem,test(1:lm))) units='US/M'; return; end
	test='RHOB';
	if(strcmp(mnem,test(1:lm))) units='K/M3'; return; end
	test='RHGF';
	if(strcmp(mnem,test(1:lm))) units='K/M3'; return; end
	test='RHGA';
	if(strcmp(mnem,test(1:lm))) units='K/M3'; return; end
	test='GRC ';
	if(strcmp(mnem,test(1:lm))) units='GAPI'; return; end
	test='GR  ';
	if(strcmp(mnem,test(1:lm))) units='GAPI'; return; end
	test='SP  ';
	if(strcmp(mnem,test(1:lm))) units='MV  '; return; end
	test='CALI';
	if(strcmp(mnem,test(1:lm))) units='MM  '; return; end
	%don't know how to do shear wave sonic
	test='DTSW';
	if(strcmp(mnem,test(1:lm))) units='US/M'; return; end
	test='NPHI';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIN';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIA';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHID';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHID';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIT';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='DPHI';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='SFLU';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='SFL ';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='ILM ';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='ILD ';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='SFLR';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='UNVI';
	if(strcmp(mnem,test(1:lm))) units='MM  '; return; end
	test='MNOR';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='MINV';
	if(strcmp(mnem,test(1:lm))) units='OHMM'; return; end
	test='DPSS';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
elseif(metric==0)
	test='DEPT';
	if(strcmp(mnem,test(1:lm))) units='FT  '; return; end
	test='AU  ';
	if(strcmp(mnem,test(1:lm))) units='US/F'; return; end
	test='DT  ';
	if(strcmp(mnem,test(1:lm))) units='US/F'; return; end
	test='SON ';
	if(strcmp(mnem,test(1:lm))) units='US/F'; return; end
	test='RHOB';
	if(strcmp(mnem,test(1:lm))) units='LB/F'; return; end
	test='RHGF';
	if(strcmp(mnem,test(1:lm))) units='LB/F'; return; end
	test='RHGA';
	if(strcmp(mnem,test(1:lm))) units='LB/F'; return; end
	test='GRC ';
	if(strcmp(mnem,test(1:lm))) units='GAPI'; return; end
	test='GR  ';
	if(strcmp(mnem,test(1:lm))) units='GAPI'; return; end
	test='SP  ';
	if(strcmp(mnem,test(1:lm))) units='MV  '; return; end
	test='CALI';
	if(strcmp(mnem,test(1:lm))) units='IN  '; return; end
	%don't know how to do shear wave sonic
	test='DTSW';
	if(strcmp(mnem,test(1:lm))) units='US/F'; return; end
	test='NPHI';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIN';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIA';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHID';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHID';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='PHIT';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='DPHI';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
	test='SFLU';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='SFL ';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='ILM ';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='ILD ';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='SFLR';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='UNVI';
	if(strcmp(mnem,test(1:lm))) units='IN  '; return; end
	test='MNOR';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='MINV';
	if(strcmp(mnem,test(1:lm))) units='OHMF'; return; end
	test='DPSS';
	if(strcmp(mnem,test(1:lm))) units='PU  '; return; end
end
