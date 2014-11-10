function prettyid=prettywid(uglyid)
% prettyid=prettywid(uglyid)
%
% PRETTYWID converts an ugly well id into a pretty one
% The prettyid will be: lsd-sect-twp-rng
% Code adapted from 
% void short_well_name_from_uwi(char* short_well_name,const char* uwi,
% 	const int form)
% found antares in /chap/extracts/prod/src/extinterp/extinterp.pc
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
if(length(uglyid)<12)
		prettyid=uglyid;
		if(length(prettyid)==0) prettyid=''; end
		return;
	end
prettyid= blanks(14);
	prettyid(1:2)=uglyid(4:5); %lsd (what a trip)
	prettyid(3)='-';
	prettyid(4:5)=uglyid(6:7); %section
	prettyid(6)='-';
	prettyid(7:9)=uglyid(8:10); %township
	prettyid(10)='-';
	prettyid(11:12)=uglyid(11:12); %range
	
