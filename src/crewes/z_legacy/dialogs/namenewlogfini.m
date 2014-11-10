function [logname,wellname,inlinecoord,zstart,zend,topsflag,logtype]=...
							namenewlogfini
% [logname,wellname,inlinecoord, zstart, zend, topsflag logtype]=namenewlogfini
%
% NAMENEWLOGFINI is called after the NAMENEWLOG dialog completes to retrieve
% the dialog responses. These are:
% logname ... string containing the name for the new log
% wellname ... string containg the name of the well or, if it is to be attached
%		to an existing well, then the number of that well
% inlinecoord ... the inline coordinate of the well. Should be ignored if it
%		is an existing well
% zstart ... starting depth to be kept for the log
% zend ... ending depth to be kept for the log
% topsflag ... if 1, then tops exist and are to be imported
%					if 0, then tops exist but are to be ignored
%					if -1, then tops don't exist
% logtype ... integer flag indicating the logtype with the following meaning
%		-1 ... unknown or un-specified
%		0  ... p-wave sonic
%		1  ... neutron density
%		2  ... formation denisty
%		3  ... gamma ray
%		4  ... spontaneous potential
%		5  ... caliper
%		6  ... s-wave sonic
%
% G.F. Margrave June 1994
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
as=get(gca,'userdata');
set(gca,'userdata',[]);
if(as==-1)
	logname=-1;
	return;
end
ind=find(isnan(as));
logname=setstr(as(1:ind(1)-1));
wellname=as(ind(1)+1:ind(2)-1);
if(length(wellname)>1)
	wellname=setstr(wellname);
end
inlinecoord=as(ind(2)+1);
zstart=as(ind(2)+2);
zend=as(ind(2)+3);
topsflag=as(ind(2)+4);
logtype=as(ind(2)+5);
