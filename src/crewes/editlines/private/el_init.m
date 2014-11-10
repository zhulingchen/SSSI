function el_init

	anchors=get(gca,'userdata');
	if(isempty(anchors)) % search the figure for the handles of all lines
		h=get(gca,'children');
		for k=1:length(h)
			if( strcmp(get(h(k),'type'),'line') )
				anchors = [anchors h(k) 0];
			end
		end
	end
	set(gca,'userdata',[]);
	% make an invisible storage bucket 
	hstore=uicontrol('style','text','visible','off','string','yrag');
	hparams=uicontrol('style','text','visible','off','string','yrag_params');
	hundo=uicontrol('style','text','visible','off','string','yrag_undo');
	hhor=uicontrol('style','text','visible','off','string','yrag_hors');
	% parameters
	% pdat(1) == xonly ... if 1 then y is not changed
	% pdat(2) == yonly ... if 1 then x is not changed
	% pdat(3) == linkmode ... if 1 thenwe are in link mode
	% pdat(4) == locate ... if 1 then we write the cursor location out
	% pdat(5) == dragmode ... if 0 then group drag is elastic, else it is constant
	% pdat(6) == smoothmode ... if 1 then we are in smooth mode
	% pdat(7) == fastopt ... if 1, then dragging a single point is done with a
	% 				less annoying display which, though less accurate, is much faster
	%				graphically
	% pdat(8) == multiple nans ... if 0 then we allow multiple nans, otherwise,
	%				we cull them
	% pdat(9) == nodelete ... if 0 then points can be delete if 1 they cannot
	%				defaults to 0
	% pdat(10) == noadd ... if 0 then points can be added, if 1 they cannot
	%				defaults to 0
	pdat=zeros(1,10);
	set(hparams,'userdata',pdat);
	
	% search for a negative anchor which indicates we should go ahead an select that
	% line for editing. Also in this loop, we create a vector of the handles alone
	% to facilitate quickly searching all lines in link mode
	% We also check to see that any handles found are children of the current
	% axes. If not, they are deleted.
	h=anchors(1);
	nh=1;
	hstart=0;
	hors=[];
	hkids=get(gca,'children');
	while( h ~= 0 )
		test=find(abs(h)==hkids);
		if(isempty(test))
			%its bogus, toss it
			npts=anchors(nh+1);
			anchors(nh:nh+1+2*npts)=[];
			if(nh>length(anchors) ) h=0;
			else h=anchors(nh);
			end
		else
			hors=[hors abs(h)];
			if( h<0 )
				if(hstart==0); % only the first negative handle is honored
					hstart=abs(h);
				end
				anchors(nh)=abs(h);
			end
			npts=anchors(nh+1);
			nh=nh+2*npts+2;
			if( nh >length(anchors) )h=0;
			else h=anchors(nh);
			end
		end
	end
	% make another storage bucket for the anchors	
	hanchors=uicontrol('style','text','visible','off','string','yrag_anchors',...
		'userdata',anchors);
		
	% put the horizon vector in its bucket
	set(hhor,'userdata',hors);

	%make sure pointer is an arrow
	set(gcf,'pointer','arrow');
		
	% start editing if requested
	if(hstart)
		set(gcf,'currentobject',hstart);
		editlines('buttondown');
	end
