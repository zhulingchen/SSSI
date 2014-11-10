function textmove(action)
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
if(nargin<1)
	set(gcf,'windowbuttondownfcn','textmove(''init'')');
	return;
end
if(strcmp(action,'init'))
	hobj=gco;
	if( strcmp(get(hobj,'type'),'text') )
			set(gcf,'windowbuttonmotionfcn','textmove(''move'')');
			set(gcf,'windowbuttonupfcn','textmove(''fini'')');
			dat1=get(hobj,'units');
			dat2=get(gcf,'units');
			set(hobj,'units','pixels');
			set(gcf,'units','pixels');
			pos=get(hobj,'position');
			xobj=pos(1);
			yobj=pos(2);
			pt=get(gcf,'currentpoint');
			xmouse=pt(1,1);
			ymouse=pt(1,2);
			dx=xobj-xmouse;
			dy=yobj-ymouse;
			set(gca,'userdata',[hobj dx dy abs(dat1)...
			abs(nan) abs(dat2)]);
			return;
	end
	return;
end
if(strcmp(action,'move'))
 objdat=get(gca,'userdata');
 hobj=objdat(1);
 dx=objdat(2);
 dy=objdat(3);
 pt=get(gcf,'currentpoint');
 pos=get(hobj,'position');
 pos(1)=pt(1,1)+dx;
 pos(2)=pt(1,2)+dy;
 set(hobj,'position',pos);
 return;
end
if(strcmp(action,'fini'))
	objdat=get(gca,'userdata');
	set(gca,'userdata',[]);
	hobj=objdat(1);
	%dat=objdat(4:length(objdat));
	%ii=find(isnan(dat));
	set(hobj,'units','pixels');
	set(gcf,'units','pixels');
	set(gcf,'windowbuttonmotionfcn','');
	set(gcf,'windowbuttonupfcn','');
	return;
end
