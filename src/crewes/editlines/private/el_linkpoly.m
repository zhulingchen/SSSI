function el_linkpoly

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
	
	% get the axes userdata which has the information on the last point
	% moved
	stuff=get(gca,'userdata');
	
	% find the last moved point
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
	hline1=dat(1);
	x1=get(hline1,'xdata');
	y1=get(hline1,'ydata');
	npts1=length(x1);
	ispoly=0;
	if(x1(1)==x1(npts1) & y1(1)==y1(npts1) )
		ispoly=1;
	end
	
	% get the last moved point
	% test it for validity
	if( it>length(x1) )
		error(' please move something on primary curve, then link');
	end
	xit=x1(it);
	yit=y1(it);

	if(~ispoly) % handle the polygon link
		npts1=npts1+1;
		x1=[x1 x1(1)];
		y1=[y1 y1(1)];
		set(dat(1),'xdata',x1,'ydata',y1);

	else % handle polygon breaking
		npts1=npts1-1;
		x1=x1(1:npts1);
		x1=[x1(it:npts1) x1(1:it-1)]; % break at 'it'
		y1=y1(1:npts1);
		y1=[y1(it:npts1) y1(1:it-1)]; % break at 'it'
		set(dat(1),'xdata',x1,'ydata',y1);
	end
