function xoff(hax)
% XOFF: turn off the x axis labels and tick marks
%
% xoff(hax)
%
% hax... axis handle. defaults to gca
if(nargin<1) hax=gca; end

set(hax,'xticklabelmode','manual','xtick',[])
