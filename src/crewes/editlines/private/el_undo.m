function el_undo

% the previous action is undone by restroing things to the state saved at the last 
% buttondown
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
				end
				if( found== 4)
					break;
				end
			end
		end
	end

	% get the undo information
	undostuff=get(hundo,'userdata');
	if(isempty(undostuff))
		return;
 end
	flag=isinf(undostuff);
	flag=find(flag==1);
	if(isempty(flag))
		otherstuff=[];
	else
		otherstuff=undostuff(flag+1:length(undostuff));
		undostuff=undostuff(1:flag-1);
	end

	% get the line data
	npts=undostuff(1);
	x=undostuff(2:npts+1);
	y=undostuff(npts+2:2*npts+1);
	ispoly=0;
	if( x(1)==x(npts) & y(1)==y(npts) )ispoly=1;end

	% get dat
	dat=undostuff(2*npts+2:length(undostuff));
	% get baddat
	baddat=get(hstor,'userdata');

	% dat is the user data of hstor and is:
	% dat(1) = handle of the selected line
	% dat(2) = original linestyle
	% dat(3) = original linewidth
	% dat(4:6) = original color
	% dat(7) = handle of duplicate line
	% dat(8) = handle of the anchor line
	% dat(9:length(dat)) = (x,y)'s of the anchors of the selected line
	if(dat(8))
		line_anchors=dat(9:length(dat));
		nanchors=length(line_anchors)/2;
	end

	set(dat(1),'xdata',x,'ydata',y);

	if( ispoly )
		set(dat(7),'xdata',x(2:npts),'ydata',y(2:npts));
	else
		set(dat(7),'xdata',x,'ydata',y);
	end
	
	if( dat(8) )
		%make sure dat(8) still exists
		hkids=get(gca,'children');
		ind=find(dat(8)==hkids);
		if(isempty(ind))
			hline3 = line(line_anchors(1:2:2*nanchors),line_anchors...
				(2:2:2*nanchors),'color',[1 0 0],'marker','o','markersize',12,...
			'linestyle','none','erasemode','xor');
			dat(8)=hline3;
		end
		set(dat(8),'xdata',line_anchors(1:2:2*nanchors),'ydata',...
			line_anchors(2:2:2*nanchors));
	elseif(baddat(8))
		delete(baddat(8));
	end

	set(hstor,'userdata',dat);

	% deal withthe otherstuff
	if(~isempty(otherstuff))
		hline=otherstuff(1);
		n2=otherstuff(2);
		x2=otherstuff(3:2+n2);
		y2=otherstuff(3+n2:2*n2+2);
		set(hline,'xdata',x2,'ydata',y2);
	end
