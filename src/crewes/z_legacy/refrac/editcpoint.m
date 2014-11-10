function points = editcpoint(action)
if(nargin<1) %set the button down function
		set(gcf,'windowbuttondownfcn','editcpoint(''init'')');
		return;
end
if(strcmp(action,'init'))
    if(strcmp(get(gcf,'selectiontype'),'normal'))
	pt=get(gca,'currentpoint');
        fprintf(1,'Creating a point');
	line(pt(1,1),pt(1,2),'color','r','linestyle','*');
	return;
    elseif(strcmp(get(gcf,'selectiontype'),'extend'))
	hpoint=gco;
	get(hpoint,'color')
	if(strcmp(get(hpoint,'type'),'line') & get(hpoint,'color') == [1 0 0])
		pt=get(gca,'currentpoint');
		set(hpoint,'userdata',pt(1,1:2));
		set(hpoint,'erasemode','xor','linestyle','*');
		
		set(gcf,'windowbuttonmotionfcn','editcpoint(''move'')');
		set(gcf,'windowbuttonupfcn','editcpoint(''fini'')');
		return;
	else
		return;
	end
    elseif(strcmp(get(gcf,'selectiontype'),'alt'))
	hpoint=gco;
	if(strcmp(get(hpoint,'type'),'line')& get(hpoint,'color') == [1 0 0])
		delete(gco);
		return;
	else
		return;
	end
    else
	return;
    end
end
if(strcmp(action,'move'))
	hpoint=gco;
	pt1=get(hpoint,'userdata');
	pt2=get(gca,'currentpoint');
	pt2=pt2(1,1:2);
		
	del=pt2-pt1;
		
	x=get(hpoint,'xdata');
	y=get(hpoint,'ydata');
		
	set(hpoint,'xdata',x+del(1));
	set(hpoint,'ydata',y+del(2));
	set(hpoint,'userdata',pt2);
		
	return;
end
	
if(strcmp(action,'fini'))
	hpoint=gco;
	set(hpoint,'erasemode','normal','linestyle','*');
		
	set(gcf,'windowbuttonmotionfcn','');
	set(gcf,'windowbuttonupfcn','');
	return;
end
if(strcmp(action,'save'))
	objs = get(gca,'children');
	n = 1;   % Counter for the number of picked points
	[number x] = size(objs);
	for i=1:number    % i is the counter for the number of axes children
		if(strcmp(get(objs(i),'linestyle'),'*'))
			point(n) = objs(i);
			get(point(n),'xdata')   % This just prints X
			n = n+1;
		end
	end
	points = point;
end
		
		
