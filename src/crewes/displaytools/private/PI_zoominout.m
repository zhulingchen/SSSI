function PI_zoominout();
% when MB1 is pressed while zoomscolling is on, user can move plot in and out
set(gcf,'windowbuttonmotionfcn','plotimage(''zoominoutmotion'')');
