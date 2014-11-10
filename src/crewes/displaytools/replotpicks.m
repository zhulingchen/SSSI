function replotpicks(figno)

if(nargin==1)
	figure(figno)
end

global PICKS PICKCOLOR

if(isempty(PICKCOLOR))
	clr='r';
else
	clr=PICKCOLOR;
end

npicks=floor(size(PICKS,1)/2);

for k=1:npicks
	line([PICKS(2*k-1,1) PICKS(2*k,1)],[PICKS(2*k-1,2) PICKS(2*k,2)],...
		[1 1],'color',clr,'linewidth',2);
end
