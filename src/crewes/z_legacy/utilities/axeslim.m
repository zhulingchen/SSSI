function axeslim(arg,pct,transfer)
% axeslim(hax,pct,transfer)
% axeslim(hax,pct)
% axeslim(hax)
%
% AXESLIM pops up a box asking for settings of the limits of the current axes.
% It the sets them when dismissed. in its simplest form, it takes a single
% argument which is a handle to the current axes:
%
% axeslim(gca)
%
% If a second argument is provided, then it gives a percentage increase for 
% the axes limits as provided by the user. In other words, if pct is 10, then
% the actual axes limits will be 10% larger than what the user actually
% provided.  This defaults to zero. 
%
% The third argument is a string with a syntactically correct Matlab command
% to be executed after the axes have been reset.
%
% G.F. Margrave, December 1993
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
if( isstr(arg) )
	action=arg;
else
	action='init';
	hax=arg;
	if(nargin < 2)
		pct=0.;
	end
	if(nargin < 3)
		transfer=[];
	end
end
if(strcmp(action,'init') )
	hfig=figure('visible','off');
	pos=get(hfig,'position');
	
	figheight=125;
	figwidth=160;
	set(hfig,'position',[pos(1),pos(2),figwidth,figheight],'visible','on');
	sep=1;
	xnow=sep;
	height=20;
	ynow=figheight-height;
	width=figwidth;
	
	hmsg=uicontrol('style','text','string','Modify the axes settings:',...
		'position',[xnow,ynow,width,height]);
	% xmin
	ynow=ynow-height;
	width=80;
	hxminlabel = uicontrol('style','text','String','Minimum x',...	
		'position',[xnow,ynow,width,height]);
	xlim=get(hax,'xlim');
	xnow=xnow+sep+width;
	hxmin = uicontrol('style','edit','String',num2str(xlim(1)),...
		'position',[xnow,ynow,width,height],'userdata',transfer);
	% xmax
	xnow=sep;
	ynow=ynow-height;
	width=80;
	hxmaxlabel = uicontrol('style','text','String','Maximum x',...	
		'position',[xnow,ynow,width,height]);
	xnow=xnow+sep+width;
	hxmax = uicontrol('style','edit','String',num2str(xlim(2)),...
		'position',[xnow,ynow,width,height]);
	% ymin
	xnow=sep;
	ynow=ynow-height;
	width=80;
	hyminlabel = uicontrol('style','text','String','Minimum y',...	
		'position',[xnow,ynow,width,height]);
	ylim=get(hax,'ylim');
	xnow=xnow+sep+width;
	hymin = uicontrol('style','edit','String',num2str(ylim(1)),...
		'position',[xnow,ynow,width,height]);
	% ymax
	xnow=sep;
	ynow=ynow-height;
	width=80;
	hymaxlabel = uicontrol('style','text','String','Maximum y',...	
		'position',[xnow,ynow,width,height]);
	xnow=xnow+sep+width;
	hymax = uicontrol('style','edit','String',num2str(ylim(2)),...
		'position',[xnow,ynow,width,height]);
	
	% done button
	xnow=sep;
	ynow=ynow-height;
	factor=1+pct/100.;
	hdone=uicontrol('style','pushbutton','string','Done','callback',...
		'axeslim(''done'')','position',[xnow,ynow,width,height],...
		'userdata',factor);
	xnow=xnow+sep+width;
	hcancel=uicontrol('style','pushbutton','string','Cancel','callback',...
		'axeslim(''cancel'')','position',[xnow,ynow,width,height]);
	% set the figures userdata
	set(gcf,'userdata',[hax hxmin hxmax hymin hymax]);
	return;
end
if( strcmp(action,'cancel') )
	close(gcf);
	return;
end
if( strcmp(action,'done') )
	h=get(gcf,'userdata');
	hdone=gco;
	factor=get(hdone,'userdata');
	hax=h(1);
	hxmin=h(2);
	hxmax=h(3);
	hymin=h(4);
	hymax=h(5);
	xmin=get(hxmin,'string');
	xmin=str2num(xmin);
	xmax=get(hxmax,'string');
	xmax=str2num(xmax);
        ymin=get(hymin,'string');
        ymin=str2num(ymin);
        ymax=get(hymax,'string');
        ymax=str2num(ymax);
	%get the old limis
	xlim=get(hax,'xlim');
	ylim=get(hax,'ylim');
	dx=(factor-1)*(xmax-xmin)/2;
	
	test=xmax-xlim(2);
	if(test)
		xmax=xmax+dx;
	end
		
	test=xmin-xlim(1);
	if(test)
		xmin=xmin-dx;
	end
	dy=(factor-1)*(ymax-ymin)/2;
	test=ymax-ylim(2);
	if(test)
		ymax=ymax+dy;
	end
	test=ymin-ylim(1);
	if(test)
		ymin=ymin-dy;
	end
	set(hax,'xlim',[xmin xmax],'ylim',[ymin ymax]);
	transfer=get(hxmin,'userdata');
	close(gcf);
	if(~isempty(transfer))
		eval(transfer);
	end
	return;
end
