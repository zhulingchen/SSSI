function savething(action)
% G.F. Margrave Feb 1994
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
if(strcmp(action,'init'))
		hax=gca;
		qs=get(hax,'userdata');
		%build the dialog box
		hdial=figure('visible','off','menubar','none');
		pos=get(hdial,'position');
		sep=1;
		nrows=5;
		%
		% assume 6 chars in 50 pixesl
		%
		charpix=6;
		%
		% now the first question
		%
		q=qs(1,:);
		ind=find(q==1);
		if(~isempty(ind))
			q=q(1:ind(1)-1);
		end
		width=50*ceil(length(q)/charpix);
		height=20;
		figheight=(height+sep)*(nrows);
		maxwidth=width;
		ynow=figheight-height-sep;
		xnow=sep;
		hq=uicontrol('style','text','string',char(q),...
					'position',[xnow ynow width height],...
					'foregroundcolor','r');
		% the yes/no toggles
		width=70;
		ynow=ynow-sep-height;
		hyes=uicontrol('style','checkbox','string','Yes','position',...
			[xnow,ynow,width,height],'callback','savething(''toggle'')');
		xnow=xnow+sep+width;
		hno=uicontrol('style','checkbox','string','No','position',...
			[xnow,ynow,width,height],'callback','savething(''toggle'')');
		if(xnow+width>maxwidth) maxwidth=xnow+width; end
		% a message
		ynow=ynow-sep-height;
		xnow=sep;
		width=200;
		hmsg=uicontrol('style','text','string','Change name if desired:',...
			'position',[xnow,ynow,width,height],'visible','off',...
			'foregroundcolor','r');
		if(xnow+width>maxwidth) maxwidth=xnow+width; end
		ynow=ynow-height-sep;
		q=qs(2,:);
		ind=find(q==1);
		if(~isempty(ind))
			q=q(1:ind(1)-1);
		end
		hname=uicontrol('style','edit','string',char(q),'position',...
			[xnow,ynow,width,height],'visible','off',...
			'backgroundcolor','c');
		%set user data on the toggles
		set(hyes,'userdata',[hno hmsg hname]);
		set(hno,'userdata',[hyes hmsg hname]);
		% ok button
		ynow=ynow-height-sep;
		width=40;
		hok=uicontrol('style','pushbutton','string','OK','position',...
			[xnow,ynow,width,height],'callback','savething(''ok'')',...);
			'foregroundcolor','r');
		% cancel button
		xnow=xnow+width+sep;
		width=60;
		hcancel=uicontrol('style','pushbutton','string','Cancel','position',...
			[xnow,ynow,width,height],'callback','savething(''ok'')');
		transfer=qs(3,:);
		ind=find(transfer==1);
		if(~isempty(ind))
			transfer=transfer(1:ind(1)-1);
		end
		set(gcf,'userdata',[hax hyes,hno,hname,hq hmsg nan abs(transfer)]);
		%get the position of the calling figure
		hparent=get(hax,'parent');
		pospar=get(hparent,'position');
		px=pospar(1)+pospar(3)/2;
		py=pospar(2)+pospar(4)/2;
		figwidth=maxwidth;
		pp=get(hq,'position');
		pp(3)=figwidth;
		set(hq,'position',pp);
		set(gcf,'position',[px py figwidth figheight],'visible','on');
		return;
	end
	%
	% toggle the yes no buttons
	%
	if(strcmp(action,'toggle'))
		htoggle=gco;
		flag=get(htoggle,'string');
		dat=get(htoggle,'userdata');
		if( strcmp(flag,'Yes') )
			set(dat(1),'value',0);
			set(dat(2),'visible','on');
			set(dat(3),'visible','on');
		else
			set(dat(1),'value',0);
			set(dat(2),'visible','off');
			set(dat(3),'visible','off');
		end
		return;
	end
	%
	% handle the OK button
	%
	if(strcmp(action,'ok'))
		h=get(gcf,'userdata');
		hbutton=gco;
		flag=get(hbutton,'string');
		if(strcmp(flag,'OK'))
			y=get(h(2),'value');
			n=get(h(3),'value');
			if( ~y & ~n )
				set(h(5),'string','Please CHOOSE Yes or No');
				return;
			end
			if( y )
				% get the name
				name=get(h(4),'string');
				%remove leading or trailing blanks
				iii=find(abs(name)~=32);
				ib=find(abs(name)==32);
				ib2=find(ib<iii(1));
				if(~isempty(name))
					name(ib(ib2))=[];
					iii=find(abs(name)~=32);
					ib=find(abs(name)==32);
					ib2=find(ib>iii(length(iii)));
					if(~isempty(name))
						name(ib(ib2))=[];
					end
				end
				if(isempty(name))
					set(h(6),'string','Provide a NON-BLANK name');
					return;
				end
			else
				name=0;
			end
		else
			name=-1;
		end
		set(h(1),'userdata',name);
		close(gcf);
	% call the transfer expression
		transfer=setstr(h(find(isnan(h))+1:length(h)));
		eval(transfer);
		return;
	end
