function el_tempfini
% this causes a temporary termination of editing. To fully terminate,
	% editlines('fini') must be called
	% reset and exit
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

	% delete any locate information that may be on
		stuff=get(gca,'userdata');
		if( length(stuff)>6 )
			htext=stuff(2,6);
		else
			htext=0;
		end
		
		if(htext)
			delete(htext);
			htext=0;
			stuff(2,6)=0;
			set(gca,'userdata',stuff);
		end

	dat=get(hstor,'userdata');

	% get the anchor info 
	anchors=get(hanchors,'userdata');

	if( length(dat)>0 )
		delete(dat(7));
		if(dat(8)) delete(dat(8)); end
		if(dat(2)==1)
			ls='-';
		elseif(dat(2)==2)
			ls='--';
		elseif(dat(2)==3)
			ls=':';
		elseif(dat(2)==4)
			ls='-.';
		elseif(dat(2)==5)
			ls='o';
		elseif(dat(2)==6)
			ls='+';
		elseif(dat(2)==7)
			ls='.';
		elseif(dat(2)==8)
			ls='*';
		elseif(dat(2)==9)
			ls='x';
		end
		
		if (dat(2) < 5)
		set(dat(1),'linestyle',ls,'linewidth',dat(3),...
			'color',dat(4:6),'erasemode','normal');
		elseif (dat(2) >= 5)
		set(dat(1),'marker',ls,'linewidth',dat(3),...
			'color',dat(4:6),'erasemode','normal');
		end
		
		% save the anchor information
		ind=find( anchors==dat(1) );
		line_anchors=dat(9:length(dat));
		nanchors=length(line_anchors)/2;
		front=anchors(1:ind);
		back=anchors(ind+2+2*anchors(ind+1):length(anchors));
		anchors=[front nanchors line_anchors back];
		if( anchors==0 )anchors=[]; end
		set(hanchors,'userdata',anchors);
		set(hstor,'userdata',[]);
	end

	%set(gca,'userdata',[]);
	set(hstor,'userdata',[]);
	set(hundo,'userdata',[]);
