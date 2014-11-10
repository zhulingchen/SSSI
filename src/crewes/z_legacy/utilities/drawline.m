function drawline(action)
% DRAWLINE works just like SELBOX except that it draws a line instead of a box
% There are alos DRAWLINEINIT and DRAWLINEFINI just like SELBOX. See help for
% SELBOX for a description
%
% by G.F. Margrave, November 1993
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
if(strcmp(action,'init') )
        p1=get(gca,'currentpoint');
        dat=get(gca,'userdata');
        if(length(dat)>=5) delete(dat(5)); end
        set(gca,'userdata',p1(1,1:2));
        set(gcf,'windowbuttonmotionfcn','drawline(''motion'')');
        return;
end
if(strcmp(action,'fini'))
        set(gcf,'windowbuttonmotionfcn','');
        return;
end
if( strcmp(action,'motion') )
% get the starting point from axes userdata
	h=get(gca,'userdata');
	
	p1=h(1:2);
	
	% delete any pre-existing line
	if(length(h)==5) delete(h(5)); end
	
% get the current point from the axes
	p2=get(gca,'currentpoint');
	p2=p2(1,1:2);
	
% draw the line
	 h=line([p1(1),p2(1)],[p1(2),p2(2)],'erasemode','xor');
	  
% update the info in userdata
 	set(gca,'userdata',[p1 p2 h]);
return;
end
