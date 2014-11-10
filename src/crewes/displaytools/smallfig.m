function smallfig(figno,figsize)
%SMALLFIG shrinks a figure to the original size
% smallfig(figno)
%
% figno is the handle of the figure to be returned to original size
% 
%
% SMALLFIG requires that both arguments be passed to the function and that
% figsize is in normalized units or pixels
%
%figsize=[lowerleftx, lowerlefty, width, height];
%
%If the function is being used in a button call back the figsize can be
%stored in the userdata of the handle for the button, then accessed by
%retriving it in the callback function as illustrated below
%
%uimenu('parent',menus, 'Label','Small Figure','Userdata',figsize,...
%           'Callback','smallfig(gcf,get(gcbo,''Userdata''))');
%
%
% H.J.E Lloyd,
scrsz = get(0,'ScreenSize');
figunits=get(figno,'Units');

if figsize(4)<=1
    sz=[figsize(1)*scrsz(3),figsize(2)*scrsz(4),figsize(3)*scrsz(3),figsize(4)*scrsz(4)];
else
    sz=figsize;
end

set(figno,'units','pixels','position',sz);
set(figno,'units',figunits);