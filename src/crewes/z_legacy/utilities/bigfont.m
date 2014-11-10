function bigfont(haxes,fontsize,fontwt,tickflag)
%
% bigfont(haxes,fontsize,fontwt)
%
% haxes ... handle of an axes object or other parent object 
%		which contains text objects
%   ******** default = gca ******
% fontsize ... desired fontsize expressed as a ratio of output
%		to input size
%   ******** default 2 ********
% fontwt ... 1 = normal
%			 2 = bold
%   *********  default = 2 *********
% tickflag ... 1 ... change tick labels
%			   0 ... dont change tick labels
%   ********* default = 1 *********
%
% Gary Margrave, CREWES
%
if(nargin<4)
	tickflag=1;
end
if(nargin<3)
	fontwt=2;
end
if(nargin<2)
	fontsize=2;
end
if(nargin<1)
	haxes=gca;
end
hkids=allchild(haxes);
for k=1:length(hkids)
	if(strcmp(get(hkids(k),'type'),'text'))
		fsize=get(hkids(k),'fontsize');
		set(hkids(k),'fontsize',fontsize*fsize);
		if(fontwt==1)
			set(hkids(k),'fontweight','normal');
		elseif(fontwt==2)
			set(hkids(k),'fontweight','bold');
		end
	end
end
if(tickflag==1)
	fsize=get(haxes,'fontsize');
	set(haxes,'fontsize',fontsize*fsize);
	if(fontwt==1)
			set(haxes,'fontweight','normal');
	elseif(fontwt==2)
			set(haxes,'fontweight','bold');
	end
end
		
