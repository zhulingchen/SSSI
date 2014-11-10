function el_buttonup(action)
	% turn off the window button functions
	set(gcf,'windowbuttonupfcn','');
	set(gcf,'windowbuttonmotionfcn','');

	
				
	% we are adding (or deleting) a point to the current line
	% or toggling the anchor status of a point
	% find the clicked point
	stuff=get(gca,'userdata');
	lit=stuff(1,1);
	it=stuff(1,2:lit+1);
	hstor=stuff(1,lit+2);
	xpt=stuff(1,lit+3);
	ypt=stuff(1,lit+4);
	hline=stuff(1,lit+5);
	hline2=stuff(1,lit+6);
	hline3=stuff(1,lit+7);
	if( hline3 )
		line_anchors=stuff(1,lit+8:length(stuff(1,:)));
		nanchors=length(line_anchors)/2;
	else
		nanchors=0;
		line_anchors=[];
	end
	xonly=stuff(2,1);
	yonly=stuff(2,2);
	killrd=stuff(2,3);
	ispoly=stuff(2,4);

	noadd=stuff(3,3);


	% get the lines data	
	x=get(hline,'xdata');
	y=get(hline,'ydata');
	npts=length(x);

	if( strcmp(action,'button1up') ) % add a point
		if(noadd)
			return;
		end
	% check to see if the point already exists
		dtest=abs(x(it)-xpt)+abs(y(it)-ypt);
		if( dtest<.5*killrd )
			return;
		end
	% find all those points which surround xpt

		ind=surround(x,xpt);

	% now for each pair of surrounding points compute the perpendicular
	% distance from xpt,ypt to the line segment connecting the pair

		m=(y(ind+1)-y(ind))./(x(ind+1)-x(ind));
		b=y(ind)-m.*x(ind);
		d=abs(m*xpt -ypt +b)./sqrt(m.*m+1);

	% find the minimum distance
		live = find(~isnan(d));
		if(isempty(live))
			editlines('tempfini');
			return;
		end
		
		it=find(d(live)==min(d(live)));
		
	% the next lines are there to trap an apparent MATLAB bug. Occaisionally, it
	% seems possible to click outside all lines and still have the current object
	% be a line. This leads to the current circumstance where we are inserting a
	% point which is way off the line. The intended behavior is a temporary fini
	% so we do that
		if(min(d)>2*killrd)
			editlines('tempfini');
			return;
		end

		% insert the point

		x=[x(1:ind(it)) xpt x(ind(it)+1:npts)];
		y=[y(1:ind(it)) ypt y(ind(it)+1:npts)];
		npts=length(x);

	elseif( strcmp(action,'button2up') )  % delete a point

		
		% make sure its not an anchor. To delete an anchor you must first free it 
		if( nanchors )
			ind=find( x(it)==line_anchors(1:2:2*nanchors) );
			if( length(ind)> 0)
				ohoh = find(y(it)==line_anchors(2*ind));
				if( length(ohoh) > 0)
					return;
				end
			end
		end
		
		%check for accidental deletion of a segment separator
		if( isnan(x(it)) | isnan(y(it)) )
			return;
		end
		
		% check for deletion of polygon endpoint
		done=0;
		if( it(1)==1 | it(1)==npts ) % it might be an end point
			if( ispoly ) % see if it is a polygon
				x=[x(2:npts-1) x(2)];
				y=[y(2:npts-1) y(2)];
				npts=length(x);
				done=1;
			end
		end
		
		if( ~done)	
			if( length(it)>1 )
				error('click closer to the point to delete it');
			end
			
			x=[x(1:it-1) x(it+1:npts)];
			y=[y(1:it-1) y(it+1:npts)];
			npts=length(x);
		end
		
	elseif( strcmp(action,'button3up') ) % toggle the points anchor status
		% see if its a polygon and we have selected an endpoint
		if(ispoly)
			if(it(1)==1 | it(1)==npts)
				it=1; % make sure only one of the endpoints is an anchor
			end
		end
		
		%check for accidental selection of a segment separator
		if( isnan(x(it)) | isnan(y(it)) )
			return;
		end
		
		if( length(it)>1 )
			it=it(1);
		end
		
		% see if it already is an anchor
		found=0; 
		if( nanchors )
			ind=find( x(it)==line_anchors(1:2:2*nanchors) );
			if( length(ind)> 0)
				ooohh = find(y(it)==line_anchors(2*ind));
				if( length(ooohh) > 0)
					% ok remove it from anchor status.
					% 2*ind(ooohh) points to the
					%  ycoord of the selected anchor
					nanchors=nanchors-1;
			 	     line_anchors=[line_anchors(1:2*ind(ooohh)-2)...
				            line_anchors(2*ind(ooohh)+1:length(line_anchors))];
					  % remove any line anchors within .5*killrd of this one
					  kmax=nanchors;
					  nanchors=0;
					  la=[];
					  for k=1:kmax
							xa=line_anchors(2*k-1);
							ya=line_anchors(2*k);
							d=sqrt( (xa-x(it))^2 + (ya-y(it))^2 );
							if( d> killrd )
								nanchors=nanchors+1;
								la=[la xa ya];
							end
						end
						line_anchors=la;
					found=1;
				end
			end
		end
		if( ~found )
			% make the point an anchor
			nanchors=nanchors+1;
			line_anchors=[line_anchors x(it) y(it)];
		end
		
		% store the anchor information
		% the user data of hstor is:
		% dat(1) = handle of the selected line
		% dat(2) = original linestyle
		% dat(3) = original linewidth
		% dat(4:6) = original color
		% dat(7) = handle of duplicate line
		% dat(8) = handle of the anchor line
		% dat(9:length(dat)) = (x,y)'s of the anchors of the selected line
		dat=get(hstor,'userdata');
		dat=[dat(1:8) line_anchors];
		set(hstor,'userdata',dat);
		
	end
				
	% check for occurance of multiple nans in a row or at the beginning &/or end of
	% the line

 singnan=stuff(3,2);
 if(singnan)
		[x,ikeep]=nanclear(x);
		if( length(x)<npts )
			npts=length(x);
			y=y(ikeep);
		end
	end

	set(hline,'xdata',x,'ydata',y);
	if( ispoly )
		set(hline2,'xdata',x(2:npts),'ydata',y(2:npts));
	else
		set(hline2,'xdata',x,'ydata',y);
	end
	
	%fart with the anchors display
	if( nanchors )
		if( hline3 ) set(hline3,'xdata',line_anchors(1:2:2*nanchors),'ydata',...
			line_anchors(2:2:2*nanchors));
		else
			hline3 = line(line_anchors(1:2:2*nanchors),line_anchors...
				(2:2:2*nanchors),'color',[1 0 0],'marker','o','markersize',12,...
			'linestyle','none','erasemode','xor');
			dat(8)=hline3;
			set(hstor,'userdata',dat);
		end
	elseif( ~nanchors & hline3 )
		delete(hline3);
		hline3=0;
		dat(8)=hline3;
		set(hstor,'userdata',dat);
	end
