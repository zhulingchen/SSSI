function el_stopmotion
 % see if a significant motion occured. If it did not, call the appropriate
 % button up function. This is done to decrease the sensitivity of the 
 % program to accidental small mouse motions
	% find the clicked point
    smoothmode=0;   % This probably isn't used anymore.  We'll define it here to make the compiler shut up.
	stuff=get(gca,'userdata');
	lit=stuff(1,1);
	it=stuff(1,2:lit+1);
	hstor=stuff(1,lit+2);
	xit=stuff(1,lit+3);
	yit=stuff(1,lit+4);
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
	locate=stuff(2,5);
	htext=stuff(2,6);
	dragmode=stuff(2,7);
	linkmode=stuff(2,8);
	fastopt=stuff(3,1);
	
	pt=get(gca,'currentpoint');
	xpt=pt(1,1);
	ypt=pt(1,2);		

	%delete any locate text that may be on
	if(htext)
		delete(htext)
		htext=0;
		stuff(2,6)=0;
		set(gca,'userdata',stuff);
	end

	% get the button
	flag=get(gcf,'selectiontype');
	if( strcmp(flag,'normal') )
		d=sqrt( (xpt-xit)^2 + (ypt-yit)^2 );
		if( d>.1*killrd )
			set(gcf,'windowbuttonmotionfcn','');
			set(gcf,'windowbuttonupfcn','');

			% test for fastopt
			if(fastopt)
				xy=get(hline2,'userdata');
				set(hline,'xdata',xy(1,:),'ydata',xy(2,:));
				if( ispoly )
                    npts = length(get(hline,'xdata'));
					set(hline2,'xdata',xy(1,2:npts),'ydata',xy(2,2:npts));  
				else
					set(hline2,'xdata',xy(1,:),'ydata',xy(2,:));
				end
			end
			return;
		else
			set(gcf,'windowbuttonmotionfcn','');
			set(gcf,'windowbuttonupfcn','');

			editlines('undo');
			if(~linkmode & ~smoothmode )
				editlines('button1up');
			elseif( linkmode==-1)
				editlines('linkpoly');
			elseif( smoothmode )
				editlines('smooth');
			end
		end
	elseif( strcmp(flag,'alt') )
		d=sqrt( (xpt-xit)^2 + (ypt-yit)^2 );
		if( d>.1*killrd )
			set(gcf,'windowbuttonmotionfcn','');
			set(gcf,'windowbuttonupfcn','');

			% test for fastopt
			if(fastopt)
				xy=get(hline2,'userdata');
				set(hline,'xdata',xy(1,:),'ydata',xy(2,:));
				if( ispoly )
                    npts = length(get(hline,'xdata'));
                    set(hline2,'xdata',xy(1,2:npts),'ydata',xy(2,2:npts));
				else
					set(hline2,'xdata',xy(1,:),'ydata',xy(2,:));
				end
			end
			return;
		else

			editlines('undo');
			if(~linkmode & ~smoothmode )
				editlines('button3up');
			elseif( linkmode)
				editlines('autoseg');
			elseif( smoothmode )
				editlines('smooth3');
			end
		end
	end

	set(gcf,'windowbuttonmotionfcn','');
    
