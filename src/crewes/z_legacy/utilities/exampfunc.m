function exampfunc(action)
%
% example matlab function
% Function reads in an ascii data file containing two columns of
% numbers [x y] and plots them. It enables zooming with mouse 
% button 1, and sets the button down function on the plotted curve
% such that if you click mouse button 3 on it it cycles through
% the axes color list. For zooming, use MB1 to draw a box or just
% click to unzoom.
% To run, just type exampfunc at matlab prompt. If you don't have
% a sample data file, read ~margrave/matlab/examp.dat which
% contains a p-wave sonic log.
%
% Suggested reading:
%	- MUG pp 2-128-> 2-149 (Introduction to programming)
%	- MUG pp 2-101-> 2-114 (Handle graphics)
%	- BAGUI pp 5-4 -> 5-8 (GUI programming method)
%	- BAGUI pp 2-1 -> 3-11 (uicontrols and menus)
%	- BAGUI pp 5-8 -> 5-20 (Handling the mouse)
%
% MUG = Matlab Users Guide
% BAGUI = Building a Graphical User Interface
% Also, look up each function you don't understand in the
% Matlab Reference Guide
%
%test number of input arguments
if( nargin < 1)
	action='init';
end
%initialize the window
if( strcmp(action,'init'))
	%make a new figure
	hfig=figure;
	%make an options menu
	hopt=uimenu(hfig,'label','Options');
	% read file option
	hread=uimenu(hopt,'label','Read File','callback',...
		'exampfunc(''read'')');
	%close the window
	hclose=uimenu(hopt,'label','Close','callback',...
		'exampfunc(''close'')');
	%save handles in figures userdata
	set(hfig,'userdata',[hread hclose]);
	return;
end
%read the file
if(strcmp(action,'read'))
	h=get(gcf,'userdata');
	hread=h(1);
	set(gcf,'pointer','watch');
	[file,path]=uigetfile('*','Select ASCII file to plot');
	if( isempty(file)|file==0 )
		disp('no file name given');
		set(gcf,'pointer','arrow');
	end
	eval(['load ' path file]);
	%load statement puts data in variable whose name is the same as
	%the file name (up to a .)
	%copy into x and y
	ind=find(file=='.');
	if(ind~=[])
		file=file(1:ind(1)-1); %shorten the name
	end
	eval(['x= ' file '(:,1);']);
	eval(['y= ' file '(:,2);']);
	%save the data and call plot
	set(hread,'userdata',[x y]);
	%plot
	exampfunc('plot');
	set(gcf,'pointer','arrow');
	return;
end
%plot the data
if(strcmp(action,'plot'))
	h=get(gcf,'userdata');
	hread=h(1);
	%get the data
	dat=get(hread,'userdata');
	x=dat(:,1);
	y=dat(:,2);
	hl=line(x,y);
	%set the button down function
	set(hl,'buttondownfcn','exampfunc(''changecolors'')','userdata',1);
	%install zooming with button 1
	simplezoom(1);
	return;
end
%close the window
if(strcmp(action,'close'))
	close(gcf);
 return;
end
if(strcmp(action,'changecolors'))
 %make sure button 3 was clicked
 flag=get(gcf,'selectiontype');
 if(~strcmp(flag,'alt'))
		return;
	end
	%get the lines handle
	hl=gco;
	%get the axes color list
	kols=get(gca,'colororder');
	nkols=size(kols,1);
	ikol=get(hl,'userdata');
	ikol=ikol+1;
	ikol=rem(ikol,nkols);
	if(ikol==0) ikol=nkols; end
	%set the lines color
	set(hl,'color',kols(ikol,:),'userdata',ikol);
	return;
end
