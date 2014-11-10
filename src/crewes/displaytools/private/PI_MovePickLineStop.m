function PI_MovePickLineStop();
	h=get(gcf,'userdata');
	delete(findobj(gcf,'type','line','tag','PICKMARKER'));
	delete(findobj(gcf,'type','text','tag','PICKTEXT'));
	h=get(gcf,'userdata');
	hzoompick=h(9);
	value=get(hzoompick,'value');
	switch value
	case 1
			selboxinit('plotimage(''zoom'')',1);
			set(gcf,'name','Seismic Image Plot, Simplezooming installed (Use MB1)');
	case 2
			drawlineinit('plotimage(''pick'')',1);
			set(gcf,'name','Seismic Image Plot, Picking resummed (Use MB1)');
	case 3
			drawlineinit('plotimage(''pick'')',1);
			set(gcf,'name','Seismic Image Plot, Picking new (Use MB1)');
			set(hzoompick,'userdata',[]);
	end
