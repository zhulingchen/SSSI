function selbox(action,button)
%
% selbox(action,button)
%
% SELBOX is somewhat like the built in RBBOX except that it 
% works in axes (data) coordinates rather than figure 
% coordinates (pixels). Thus the selection rectangle can 
% be used directly in data coordinates to test data for
% inclusion. SELBOX uses the userdata of the current axes
% for storage and assumes that it initially contains [pt1x pt1y] the
% x and y coordinates of the starting point of the selection rectangle. 
% So that SELBOX can be convieniently called several times in a row
% without outside intervention,SELBOX begins by deleting the graphics
% object whose handle is the fifth entry in the axes user data.
% An abort will occur if this user data entry is non-null with an
% invalid graphics handle.) After returning, the axes user data will
% contain the x,y coordinate of the 2 points defining the 
% rectangle and the graphics handle of the selection rectangle: 
% [p1x p1y p2x p2y handle]
%
% Several convienience functions make it easier to use SELBOX and
% you actually will rarely need to call SELBOX directly. SELBOXINIT
% should be called to begin the selection box process and SELBOXFINI
% to end it. You need do nothing else. SELBOXFINI, will return five
% values: [pt1x pt1y pt2x pt2y hbox]
% Where pt1 and pt2 are the points defining the rectangle and hbox
% is the graphics handle to its on-screen image. (Note that 
% SELBOXFINI needs no arguments while SELBOXINIT has one optional
% argument described below)
% Example:
% [x,y,z]=peaks; % get something to plot
% pcolor(x,y,z); % plot it
% 
% selboxinit % begin the selection box
% % the selection box is drawn on screen
% vals=selboxfini; % obtain the selection box coordinates
%
% Often you will want to do some action of your own upon a mouse button up
% event. This is easily done by passing a string argument to SELBOXINIT
% which contains any Matlab command (or series of commands).
% selboxinit('do_my_zoom')
%
% SELBOXINIT now takes a second argument which is the button number
% to be used in drawing the zoom box. Drawing with any other button has
% no effect. Default is button 1.
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
	if(nargin<2)
		button=1;
	end
if(strcmp(action,'init') )
		% get the buttion
		go=0;
		flag=get(gcf,'selectiontype');
		if( strcmp(flag,'normal') & button==1)
			go=1;
		elseif( strcmp(flag,'extend') & button==2)
			go=1;
		elseif( strcmp(flag,'alt') & button==3)
			go=1;
		end
		if(go)
		  p1=get(gca,'currentpoint');
		  dat=get(gca,'userdata');
		  if(length(dat)>=5) delete(dat(5)); end
		  set(gca,'userdata',p1(1,1:2));
		  set(gcf,'windowbuttonmotionfcn','selbox(''motion'')');
		end
		return;
end
if(strcmp(action,'fini'))
		% get the buttion
		go=0;
		flag=get(gcf,'selectiontype');
		if( strcmp(flag,'normal') & button==1)
			go=1;
		elseif( strcmp(flag,'extend') & button==2)
			go=1;
		elseif( strcmp(flag,'alt') & button==3)
			go=1;
		end
		if(go)
			set(gcf,'windowbuttonmotionfcn','');
			%     check for valid zoom box
			zoombox=get(gca,'userdata');
			if(length(zoombox) < 5)
				xlim=get(gca,'xlim');
				ylim=get(gca,'ylim');
				p1(1)=xlim(1);p1(2)=ylim(1);
				p2(1)=xlim(2);p2(2)=ylim(2);
				h=line([p1(1),p2(1),p2(1),p1(1),p1(1)],[p1(2),p1(2),p2(2),p2(2),...
					p1(2)],'erasemode','xor','linestyle','-.','linewidth',4,...
					'color',[.5 .5 .5]);
				set(gca,'userdata',[p1 p2 h]);
			end
		end
       
return;
end
if( strcmp(action,'motion') )
% get the starting point from axes userdata
	h=get(gca,'userdata');
	
	p1=h(1:2);
	
	% delete any pre-existing rectangle
	if(length(h)==5) delete(h(5)); end
	
% get the current point from the axes
	p2=get(gca,'currentpoint');
	p2=p2(1,1:2);
% see if the box is too small
	d=sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
	if( d<100*eps )
		xlim=get(gca,'xlim');
		ylim=get(gca,'ylim');
		p1(1)=xlim(1);p1(2)=ylim(1);
		p2(1)=xlim(2);p2(2)=ylim(2);
	end
% draw the box
	 h=line([p1(1),p2(1),p2(1),p1(1),p1(1)],[p1(2),p1(2),p2(2),p2(2),...
	p1(2)],'erasemode','xor','linestyle','-.','linewidth',4,...
	'color',[.5 .5 .5]);
%h=line([p1(1),p21),p2(1),p1(1),p1(1)],[p1(2),p1(2),p2(2),p2(2),...
%	p1(2)],'erasemode','xor','linestyle','-.','linewidth',4);
% update the info in userdata
 	set(gca,'userdata',[p1 p2 h]);
	return;
end
