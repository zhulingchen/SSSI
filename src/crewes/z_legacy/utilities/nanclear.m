function [xout,ikeep]=nanclear(xin)
%
% [xout,ikeep]=nanclear(xin)
%
% NANCLEAR is used by LOGSEC to clear up un-needed NAN's in a multisegmented
% horizon. These are defined as NAN's aat the beginning or end of the 
% horizon and any occurrance of multiple NAN's in a row. Multiple NAN's 
% are reduced to a single NAN. ikeep is a vector of indicies such that
% xout=xin(ikeep)
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
id=isnan(xin);
idead=find(id);
ikeep=1:length(xin);
% simple return for no-nans
if( length(idead) == 0 )
	xout=xin;
	return;
end
% check for non-isolated nans
it=id+[0 id(1:length(id)-1)];
ikeep=find(it<2);
% toss any nans at the beginning or end
if(isnan(xin(ikeep(1))))
	ikeep=ikeep(2:length(ikeep));
end
if(isnan(xin(ikeep(length(ikeep)))))
	ikeep=ikeep(1:length(ikeep)-1);
end
xout=xin(ikeep);
