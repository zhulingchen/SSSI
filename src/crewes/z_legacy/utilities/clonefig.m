function hclone=clonefig(arg,arg2,arg3,arg4)

% hclone =clonefig(hfig,buttonnum,lwfact,titflag)
%	(all parameters may be defaulted)
%
% CLONEFIG is designed to produce a copy of a figure which contains any
% number of 2-D axes with any number of lines and text. The copy is the
% same as the original in color, style, and size of lines and text with
% the possible exception of actions provided by the lwfact parameter.
% The cloned figure provides facilities for scaled hardcopy and zooming
% regardless of what was available on the original figure. Additionally,
% if the line objects have userdata which is a string, then their button
% down function is set to display this string in a message window. This
% is a standard mechanism allowing lines to identify themselves.
%	hfig ... handle of the figure to be cloned
%			******** default is current figure *********
%	buttonnum ... mouse button number which will be used for zooming
%			******** default is 3 ********
%	lwfact ... multiplicative factor to be applied to the linewidths
%		of lines as they are reproduced. A value of 1 will cause the linewidths
%		to be the same as hfig while 2 will double all linewidths.
%			******** default is 1 *******
%	titflag ... if 1, then the title of hfig will be modified to be prepended
%		with the string 'Clone of:'. If 0, the title is unchanged.
%			******** default is 1 *******
%	hclone ... handle of the clone figure
%
% G.F. Margrave Summer 1994
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

if(nargin<1)
	action='init';
	hfig=gcf;
	if(nargin<2)
		button=3;
	else
		button=arg2;
	end
elseif( ~isstr(arg))
	hfig=arg;
	action='init';
	if(nargin<2)
		button=3;
	else
		button=arg2;
	end
	if( nargin < 3)
		lwfact=1;
	else
		lwfact=arg3;
	end
	if( nargin < 4)
		titflag=1;
	else
		titflag=arg4;
	end
else
	action=arg;
end

if(strcmp(action,'init') )
	pos=get(hfig,'position');
	pos(1)=pos(1)+.5*pos(3);
	pos(2)=pos(2)-.5*pos(4);
	hclone=figure('position',pos,'menubar','none');
	
	%get the children of the figure
	hfigkids=get(hfig,'children');
	
	%loop over and act on any that are axes
	for k=1:length(hfigkids)
		if( strcmp(get(hfigkids(k),'type'),'axes') )
			
			hax=hfigkids(k);
			%make an axes of this size in the new figure
			axpos=get(hax,'position');
			axun=get(hax,'units');
			xdir=get(hax,'xdir');
			xgrid=get(hax,'xgrid');
			xlabel=get(hax,'xlabel');
			xlim=get(hax,'xlim');
			xlimmode=get(hax,'xlimmode');
			xtick=get(hax,'xtick');
			xticklabel=get(hax,'xticklabel');
			ydir=get(hax,'ydir');
			ygrid=get(hax,'ygrid');
			ylabel=get(hax,'ylabel');
			ylim=get(hax,'ylim');
			ylimmode=get(hax,'ylimmode');
			ytick=get(hax,'ytick');
			yticklabel=get(hax,'yticklabel');
			
			hnewax=axes('position',axpos,'units',axun);
			
			set(hnewax,'xdir',xdir,'xgrid',xgrid);
			set(hnewax,'xlabel',text(0,0,get(xlabel,'string')),'xlim',xlim,'xlimmode',xlimmode);
			set(hnewax,'xtick',xtick,'xticklabel',xticklabel);
			set(hnewax,'ydir',ydir,'ygrid',ygrid);
			set(hnewax,'ylabel',text(0,0,get(ylabel,'string')),'ylim',ylim,'ylimmode',ylimmode);
			set(hnewax,'ytick',ytick,'yticklabel',yticklabel);

				
			htit=get(hax,'title');
			str=get(htit,'string');
			if(~strcmp(str,''))
				if(titflag)
					title(['Clone of: ' str]);
				else
					title(str);
				end
			else
				if(titflag)
					title('Cloned Figure');
				end
			end
			
			haxkids=get(hax,'children');
			
			%loop over axes children and copy any text and line objects
			for kk=fliplr(haxkids')
			
				if( strcmp(get(kk,'type'),'line') )
					
					x=get(kk,'xdata');
					y=get(kk,'ydata');
					kol=get(kk,'color');
					ls=get(kk,'linestyle');
					lw=lwfact*get(kk,'linewidth');
					
					ud=get(kk,'userdata');
					if( ~isstr(ud) )
						ud=[];
					end
					
					if(isempty(ud))
						line(x,y,'color',kol,'linestyle',ls,'linewidth',lw);
					else
						line(x,y,'color',kol,'linestyle',ls,'linewidth',lw,...
							'userdata',ud,'buttondownfcn','clonefig(''sayhey'')');
					end
						
				elseif( strcmp(get(kk,'type'),'text') )
					str=get(kk,'string');
					pos=get(kk,'position');
					un=get(kk,'units');
					kol=get(kk,'color');
					ha=get(kk,'horizontalalignment');
					va=get(kk,'verticalalignment');
					r=get(kk,'rotation');
					
					text('position',pos,'string',str,'units',un,'color',kol,...
					'horizontalalignment',ha,'verticalalignment',va,...
					'rotation',r);
				elseif( strcmp(get(kk,'type'),'patch') )
					x=get(kk,'xdata');
					y=get(kk,'ydata');
					fkol=get(kk,'facecolor');
					ekol=get(kk,'edgecolor');
					cdata=get(kk,'cdata');
					ud=get(kk,'userdata');

					if( ~isstr(ud) )
						ud=[];
					end
					
					if(isempty(ud))
						patch('xdata',x,'ydata',y,'edgecolor',ekol,...
							'facecolor',fkol,'cdata',cdata);
					else
						patch('xdata',x,'ydata',y,'edgecolor',ekol,...
							'facecolor',fkol,'cdata',cdata,...
							'userdata',ud,'buttondownfcn','clonefig(''sayhey'')');
					end

				end
			end
			
		end
	end
	%add a close button, a hardcopy button, a msg panel, and install simple zooming

    
    
    
    
    
	sep=1;
	xnow=sep;
	ynow=sep;
	width=40;
	height=20;
	hclose=uicontrol('style','pushbutton','string','Close','position',...
		[xnow,ynow,width,height],'callback','close(gcf)');
	
	xnow=xnow+width+sep;
	width=80;	
	hhardcopy=uicontrol('style','pushbutton','string','Hardcopy','position',...
		[xnow,ynow,width,height],'callback','clonefig(''hardcopy'')',...
		'userdata',[0 0]);
		
	xnow=xnow+width+sep;
	width=400;
	if(button==1)
		msg= 'MB1 drag -> zoom  ... MB1 click -> unzoom';
	elseif(button==2)
		msg= 'MB2 drag -> zoom  ... MB2 click -> unzoom';
	else
		msg= 'MB3 drag -> zoom  ... MB3 click -> unzoom';
	end
	hmsg=uicontrol('style','text','string',msg,...
		'position',[xnow,ynow,width,height]);
		
	set(gcf,'userdata',[hclose,hhardcopy,hmsg]);
	
	simplezoom(button);
	
	return;

end

if(strcmp(action,'hardcopy') )
	h=get(gcf,'userdata');
	hhardcopy=h(2);
	
	scales=get(hhardcopy,'userdata');

	shardcopy(gcf,'clonefig(''hardcopy2'')',scales(1),scales(2),...
		'temp.ps');

	return;
end

if(strcmp(action,'hardcopy2') )
	h=get(gcf,'userdata');
	hhardcopy=h(2);
	hmessage = h(3);
	
	% get the dialog answers and test for reasonableness
	[ps_scale,xlength,ylength,xscale,yscale]=shardcopyfini;

	if( ps_scale== -1) %test for a cancel
		set(hmessage,'string','Plot cancelled');
		return;
	end
	
	%remember the scales
	set(hhardcopy,'userdata',[xscale yscale]);

	%put out a message
	set(hmessage,'string',...
		['plotsize is ' num2str(xlength) ' by ' num2str(ylength) ...
		' inches. use ps_scale = ' num2str(ps_scale) ' in CHVSUB']);
	
	return
end
		
if(strcmp(action,'sayhey'))
	h=get(gcf,'userdata');
	hmsg=h(3);

	dat=get(gco,'userdata');

	if(isstr(dat))
		set(hmsg,'string',strunpad(dat));
	end

	return;
end
