function refresh(h)
% REFRESH REFRESH causes the current figure window to be
%	redrawn. REFRESH(h) causes the figure window specified
%	by handle h to be redrawn.
%	D. Thomas   5/26/93
%       Copyright (c) 1984-93 by The MathWorks, Inc.
if nargin==1,
	if strcmp(get(h,'type'),'figure'),
		set(h,'color',get(h,'color'))
	else,
		disp('Handle does not refer to a figure object')
	end
else,
	set(gcf,'color',get(gcf,'color'))
end
