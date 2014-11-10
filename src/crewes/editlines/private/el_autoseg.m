function el_autoseg
	% turn off the window button functions
	set(gcf,'windowbuttonupfcn','');
	set(gcf,'windowbuttonmotionfcn','');
	
	% Get the storage buckets
	h=get(gcf,'children');
	found=0;
	for k=1:length(h)
		if( strcmp(get(h(k),'type'),'uicontrol') )
			if( strcmp(get(h(k),'style'),'text') )
				if( strcmp(get(h(k),'string'),'yrag') )
					hstor=h(k);
					found=found+1;
				elseif( strcmp(get(h(k),'string'),'yrag_params') )
					hparams=h(k);
					found=found+1;
				elseif( strcmp(get(h(k),'string'),'yrag_anchors') )
					hanchors=h(k);
					found=found+1;
				elseif( strcmp(get(h(k),'string'),'yrag_undo') )
                	hundo=h(k);
                	found=found+1;
                elseif( strcmp(get(h(k),'string'),'yrag_hors') )
                	hhors=h(k);
                	found=found+1;
				end
				if( found== 5)
					break;
				end
			end
		end
	end

	dat=get(hstor,'userdata');
	% the user data of hstor is:
	% dat(1) = handle of the selected line
	% dat(2) = original linestyle
	% dat(3) = original linewidth
	% dat(4:6) = original color
	% dat(7) = handle of duplicate line
	% dat(8) = handle of the anchor line
	% dat(9:length(dat)) = (x,y)'s of the anchors of the selected line
	
	% get the axes userdata which has the information on the clicked point
	stuff=get(gca,'userdata');
	
	% find the clicked point
	stuff=get(gca,'userdata');
	lit=stuff(1,1);
	it=stuff(1,2:lit+1);
	hline3=stuff(1,lit+7);
	if( hline3 )
		line_anchors=stuff(1,lit+8:length(stuff(1,:)));
		nanchors=length(line_anchors)/2;
	else
		nanchors=0;
		line_anchors=[];
	end

	%get the data from the current curve
	hline=dat(1);
	x=get(hline,'xdata');
	y=get(hline,'ydata');
	npts=length(x);
	ispoly=0;
	if(x(1)==x(npts) & y(1)==y(npts) )
		ispoly=1;
	end
	
	% get the clicked point
	xit=x(it);
	yit=y(it);
	
	% return if the point is the first or last
	if( it==1 | it==npts )
		return;
	end
	
	% make sure it is not already segmented
	if( isnan(x(it+1)) | isnan(x(it-1)) )
		return;
	end
	
	% duplicate the point and put a nan between
	
	x=[x(1:it-1) xit nan xit x(it+1:npts)];
	y=[y(1:it-1) yit nan yit y(it+1:npts)];
	
	% see if the point is already an anchor
		
	done=0;
	if(nanchors)
		ind=find(xit==line_anchors(1:2:2*nanchors));
		if(length(ind)>0)
			ia2=find(yit==line_anchors(2*ind));
			if(length(ia2)>0)
					done=1;
			end
		end
	end
	
	if(~done)
			nanchors=nanchors+1;
			line_anchors=[line_anchors xit yit];
			% line_anchors get updated
			dat=[dat(1:8) line_anchors];
			set(hstor,'userdata',dat);
	end

	% redisplay the line

	set(hline,'xdata',x,'ydata',y);% the current line

	if( ispoly )
			set(dat(7),'xdata',x(2:npts),'ydata',y(2:npts));
	else
			set(dat(7),'xdata',x,'ydata',y);
	end
	
	
	%fart with the anchors display
	if( nanchors )
		if( dat(8) ) set(dat(8),'xdata',line_anchors(1:2:...
				2*nanchors),'ydata',line_anchors(2:2:2*nanchors));
		else
			dat(8) = line(line_anchors(1:2:2*nanchors),...
					line_anchors(2:2:2*nanchors),'color',...
					[1 0 0],'marker','o','markersize',12,...
					'linestyle','none','erasemode','xor');
			set(hstor,'userdata',dat);
		end
	end
