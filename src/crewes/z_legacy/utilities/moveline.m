function moveline(action)
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
if(nargin<1) %set the button down function
		set(gcf,'windowbuttondownfcn','moveline(''init'')');
		return;
	end
	if(strcmp(action,'init'))
		hline=gco;
		if(~strcmp(get(hline,'type'),'line'))
			return;
		end
		pt=get(gca,'currentpoint');
		set(hline,'userdata',pt(1,1:2));
		set(hline,'erasemode','xor','linestyle','.');
		
		set(gcf,'windowbuttonmotionfcn','moveline(''move'')');
		set(gcf,'windowbuttonupfcn','moveline(''fini'')');
		return;
	end
	if(strcmp(action,'move'))
		hline=gco;
		
		pt1=get(hline,'userdata');
		pt2=get(gca,'currentpoint');
		pt2=pt2(1,1:2);
		
		del=pt2-pt1;
		
		x=get(hline,'xdata');
		y=get(hline,'ydata');
		
		set(hline,'xdata',x+del(1));
		set(hline,'ydata',y+del(2));
		set(hline,'userdata',pt2);
		
		return;
	end
	
	if(strcmp(action,'fini'))
		hline=gco;
		set(hline,'erasemode','normal','linestyle','-');
		
		set(gcf,'windowbuttonmotionfcn','');
		set(gcf,'windowbuttonupfcn','');
		return;
	end
		
		
