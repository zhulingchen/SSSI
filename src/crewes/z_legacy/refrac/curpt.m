function curpt(action)
 if(nargin<1) %set the button down function
		set(gcf,'windowbuttondownfcn','curpt(''init'')');
		return;
	end
	if(strcmp(action,'init'))
		%hline=gco;
		%if(~strcmp(get(hline,'type'),'line'))
		%	return;
		%end
		pt=get(gca,'currentpoint');
		%set(hline,'userdata',pt(1,1:2));
		%set(hline,'erasemode','xor','linestyle','.');
		
		%set(gcf,'windowbuttonmotionfcn','curpt(''move'')');
		%set(gcf,'windowbuttonupfcn','curpt(''fini'')');
		line(pt(1,1),pt(1,2),'color','r','linestyle','*');
		return;
	end
	if(strcmp(action,'move'))
		hline=gco;
		
		pt1=get(hline,'userdata');
		pt2=get(gca,'currentpoint');
		pt2=pt2(1,1:2);
		
		del=pt2-pt1;
		
		x=get(hline,'xdata');
		y=get(hline,'ydata');
		
		set(hline,'xdata',x+del(1));
		set(hline,'ydata',y+del(2));
		set(hline,'userdata',pt2);
		
		return;
	end
	
	if(strcmp(action,'fini'))
		hline=gco;
		set(hline,'erasemode','normal','linestyle','-');
		
		set(gcf,'windowbuttonmotionfcn','');
		set(gcf,'windowbuttonupfcn','');
		return;
	end
end
		
		
