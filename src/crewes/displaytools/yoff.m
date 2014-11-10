function yoff(hax)
% YOFF: turn off the y axis labels and tick marks
%
% yoff(hax)
%
% hax... axis handle. defaults to gca
if(nargin<1) hax=gca; end

set(hax,'yticklabelmode','manual','ytick',[])
