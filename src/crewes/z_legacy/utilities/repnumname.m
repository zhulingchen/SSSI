function [names,newname]=repnumname(names,oldname,newname)
%
% [names,newname]=addnumname(names,oldname,newname)
%
% Given a string matrix of possibly numbered names (see NUMBERNAMES)
% repnumname finds an old numbered name and replaces it with a new name 
% which is numbered as needed
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
	ln=length(newname);
	lo=length(oldname);
	[n,c]=size(names);
	
	if(ln>c-2)
		error(' new name too long ');
	end
	
	%find the old name
	indold=strmatfind(names,oldname);
	
	
	%find any names matching the newname
	ind=strmatfind(names(:,1:c-2),newname);
	
	%replace without index
	newname=[newname blanks(c-ln)];
	names(indold,:)=newname;
	
	if(~isempty(ind))
		num=-1;
		for k=ind
			% determine index
			anum=str2num(names(k,c-1:c));
			if(isempty(anum)) anum=1; end
			if(anum>num) num=anum; end
			
		end
		num=num+1;
		snum=int2str(num);
		ls=length(snum);
		
		names(indold,c-ls+1:c)=snum;
		newname(c-ls+1:c)=snum;
	end
