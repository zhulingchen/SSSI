% Clears and/or delete all figure children of type 'uicontrol' or 'axes'
%
% action = 
function clearkids(action)
c = get(gcf, 'children');
for i=1:length(c)
  if( strcmp( get(c(i), 'type'), 'uicontrol') )
     delete(c(i));
  elseif( strcmp( get(c(i), 'type'), 'axes') )
     axes(c(i));
     cla;
  end
end
