function [shots, currentshot, axeslist] = cvpinfo(action, shots, currentshot, axeslist )
if( nargin < 1 )
   action = 'get';
end
% Get the handle of the 'set range' uicontrol
c = get(gcf, 'children');
[nc tmp] = size(c);
for i=1:nc
   if( strcmp(get(c(i),'type'),'uicontrol') )
      if( strcmp(get(c(i),'string'), 'Set range'))
         rangeb = c(i);
      elseif( strcmp(get(c(i),'string'), 'Next pair'))
         nextb = c(i);
      elseif( strcmp(get(c(i),'string'), 'Previous pair'))
         prevb = c(i);
      end
   end
end
if(isempty(rangeb))
   error('cvpinfo: did not find set range button!');
end
if(isempty(nextb))
   error('cvpinfo: did not find next pair button!');
end
if(isempty(prevb))
   error('cvpinfo: did not find previous pair button!');
end
if( strcmp(action, 'get') )
   shots = get(rangeb,'userdata');
   currentshot = get(nextb,'userdata');
   axeslist = get(prevb,'userdata');
end
if( strcmp(action, 'set') )
   if( nargin < 4 )
      error('cvpinfo: not enough input arguements');
   end
   set(rangeb, 'userdata', shots)
   set(nextb, 'userdata', currentshot);
   set(prevb, 'userdata', axeslist );
end
