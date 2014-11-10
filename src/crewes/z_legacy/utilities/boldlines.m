function boldlines(haxes,linewidth,markersize)
%
% boldlines(haxes,linewidth,markersize)
%
% haxes ... handle of an axes object or other parent object 
%		which contains line objects
%   ******** default = gca ******
% linewidth ... desired line width expressed as a ratio of output
%		to input size
%   ******** default 4 ********
% markersize ... desired line width expressed as a ratio of output
%		to input size
%   ******** default 2 ********
%
% Gary Margrave, CREWES
%

if(nargin<3)
	markersize=2;
end
if(nargin<2)
	linewidth=4;
end
if(nargin<1)
	haxes=gca;
end

hkids=allchild(haxes);

for k=1:length(hkids)
	if(strcmp(get(hkids(k),'type'),'line'))
		lsize=get(hkids(k),'linewidth');
		set(hkids(k),'linewidth',linewidth*lsize);
		msize=get(hkids(k),'markersize');
		set(hkids(k),'markersize',markersize*msize);
	end
end
