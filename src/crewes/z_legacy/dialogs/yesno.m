function yesno(action)
%
%
% G.F. Margrave Feb 1994
%

if(strcmp(action,'init'))
		hax=gca;
		qs=get(hax,'userdata');

		%build the dialog box
		hdial=figure('visible','off','menubar','none');
		pos=get(hdial,'position');
		sep=1;
		nrows=2;

		%
		% assume 6 chars in 50 pixesl
		%
		charpix=6;
		%
		% now the question
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
					'position',[xnow ynow width height]);

		% the yes/no buttons

		width=40;
		ynow=ynow-sep-height;

		hyes=uicontrol('style','pushbutton','string','Yes','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');
		xnow=xnow+sep+width;

		hno=uicontrol('style','pushbutton','string','No','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');

		
		% cancel button

		xnow=xnow+sep+width;
		width=60;
		hcancel=uicontrol('style','pushbutton','string','Cancel','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');

		if(xnow+width>maxwidth) maxwidth=xnow+width; end

		transfer=qs(3,:);
		ind=find(transfer==1);
		if(~isempty(ind))
			transfer=transfer(1:ind(1)-1);
		end
		set(gcf,'userdata',[hax nan abs(transfer)]);

		%get the position of the calling figure
		hparent=get(hax,'parent');
		pospar=get(hparent,'position');
		figwidth=maxwidth;
		px=pospar(1)+pospar(3)/2-figwidth/2;
		py=pospar(2)+pospar(4)/2-figheight/2;
		

		set(gcf,'position',[px py figwidth figheight],'visible','on');

		return;

	end

	%
	% handle a button
	%
	if(strcmp(action,'button'))
		h=get(gcf,'userdata');

		hbutton=gco;

		flag=get(hbutton,'string');

		if(strcmp(flag,'Yes'))
			reply=1;
		elseif(strcmp(flag,'No'))
			reply=0;
		else
			reply=-1;
		end
		
		set(h(1),'userdata',reply);

		close(gcf);

	% call the transfer expression
		transfer=setstr(h(find(isnan(h))+1:length(h)));

		eval(transfer);

		return;

	end
