function simplezoom(arg,arg2)
% simplezoom(button,transfer)
%
% Installs simple zooming action on the current figure. To
% Zoom any 2-D axes, simply make sure the axes parent figure is
% current and type: simplezoom
% at the matlab prompt to install the zooming function. Then simply
% draw a rectangle by clicking and dragging the mouse. The axis will
% 'zoom' at the completion of the draw. To unzoom, enter a rectangle
% of zero size which is done by a simple mouse click without any 
% motion.
% button ... refers to the mouse button to be used to draw
%		the zoom box. Any other button will cause no action.
%  ********* Button default is 1 ***********
% Transfer ... is a string containing any legal matlab command
% 		which is to be called at the completion of the zoom. 
%		This is key to using SIMPLEZOOM in a larger program so that
%		control can be returned to the larger program following
%		the zoom.
% ********** Default is '' ************
%
% G.F. Margrave, Feb 1994
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
	arg2='';
else
	button=arg2;
end
if(nargin<1)
	arg=1;
else
	action=arg;
end
if( ~isstr(arg) )
	% initialize
	%
	if(~isstr(arg2)) cmd=''; else cmd=arg2; end
	if(arg==1)
		selboxinit(['simplezoom(''zoomit'',1);' arg2],arg);
	elseif(arg==2)
		selboxinit(['simplezoom(''zoomit'',2);' arg2],arg);
	elseif(arg==3)
		selboxinit(['simplezoom(''zoomit'',3);' arg2],arg);
	end
	return;
end
if( strcmp(action,'zoomit') )
 % determine the button type
 flag=get(gcf,'selectiontype');
 go=0;
 if( strcmp(flag,'normal') & button==1)
		go=1;
	elseif( strcmp(flag,'extend') & button==2)
		go=1;
	elseif( strcmp(flag,'alt') & button==3)
		go=1;
	end
	if(go)
		box=selboxfini;
		if(length(box)==5)
			delete(box(5));
		end
		xmin=min([box(1) box(3)]);
		xmax=max([box(1) box(3)]);
		ymin=min([box(2) box(4)]);
		ymax=max([box(2) box(4)]);
		%get the current axis settings
		xlim=get(gca,'xlim');
		ylim=get(gca,'ylim');
		test1=xmin-xlim(1)+xmax-xlim(2)+ymin-ylim(1)+ymax-ylim(2);
		test2=(xmin-xmax)*(ymin-ymax);
		if(abs(test1)<10*eps | abs(test2)< 10*eps)
			axis('auto')
		else
			set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax]);
		end
	end
	return;
end
