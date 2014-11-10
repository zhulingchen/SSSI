function bigfig(figno)
%BIGFIG enlarges a figure to an optimal size for a slide
% bigfig(figno)
%
% figno is the handle of the figure to be enlarged
% figno defaults to gcf
%
% Normally, BIGFIG enlarges the figure to the size of your screen
% However, you can control this behavior by defining some globals
% BIGFIG_X BIGFIG_Y BIGFIG_WIDTH and BIGFIG_HEIGHT. In order, these are 
% the x and y coordiates of the lower left corner and the width and 
% height of the figure after enlargement. Specify them in pixels.
% They default to 1,1,screen width and screen height. A good place 
% to define them is in your startup.m file
%
% G.F. Margrave, CREWES
%
if(nargin<1); figno=gcf; end
global BIGFIG_X BIGFIG_Y BIGFIG_WIDTH BIGFIG_HEIGHT
scr=get(0,'screensize');
units=get(figno,'Units');
if(isempty(BIGFIG_X))
    x=scr(1);
else
    x=BIGFIG_X;
end
if(isempty(BIGFIG_Y))
    y=scr(2);
else
    y=BIGFIG_Y;
end
if(isempty(BIGFIG_WIDTH))
    w=scr(3);
else
    w=BIGFIG_WIDTH;
end
if(isempty(BIGFIG_HEIGHT))
    h=scr(4);
else
    h=BIGFIG_HEIGHT;
end
    
set(figno,'Units','Pixels','position',[x y w h]);
set(figno,'Units',units);