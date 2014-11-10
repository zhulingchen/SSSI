function lash=lashtops(lash,topnames,ztops)
% lash=lashtops(lash,topnames,ztops)
%
% LASHTOPS searches an LAS header for a tops block and replaces it
% with one containing the topnames and tops given. It also searches
% for and eliminates multiple tops blocks (as were occaisionally
% created by logedit.)
%
% lash ... las header as a string matrix
% tops ... string matrix of tops names
% ztops ... vector of tops (numerical formation depths)
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
%first eliminate any tops blocks
[tblk,l1,l2]=lashgetblk(lash,'tops');
while(~isempty(tblk))
	lash(l1:l2,:)=[];
	[tblk,l1,l2]=lashgetblk(lash,'tops');
end
if(~isempty(ztops))
	%make new tops block
	n=length(ztops);
	stops=32*ones(n,20);
	lm=1;
	for k=1:length(ztops)
		s=sprintf('%g',ztops(k));
		l=length(s);
		stops(k,1:l)=s;
		if(l>lm) lm=l; end
	end
	stops=setstr(stops(:,1:lm));
	
	%make comment line
	ln=size(topnames,2);
	cmnt='#TOPS NAME          .        DEPTH:';
	if(ln<20)
		topnames=[topnames setstr(32*ones(n,20-ln-1))];
	end
		
	tblk=lashblkcreate('tops',cmnt,2,topnames,[],stops,[]);
	
	%insert
	lash=lashsetblk(lash,'tops',tblk);
end
	
	
