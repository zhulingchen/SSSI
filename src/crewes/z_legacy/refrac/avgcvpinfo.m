function [shotrange, shotlist, currentshot, axeslist, stddev] = ...
         avgcvpinfo(action, shotrange, shotlist, currentshot, axeslist, stddev)
% shotrange   - 
% shotlist    - 
% currentshot - 
% axeslist    - 
% stddev      -
if( nargin < 1 )
   action = 'get';
end
% Get the handle of the 'set range' uicontrol
c = get(gcf, 'children');
[nc tmp] = size(c);
for i=1:nc
   if( strcmp(get(c(i),'type'),'uicontrol') )
      if( strcmp(get(c(i),'string'), 'Set shot range'))
         rangeb = c(i);
      elseif( strcmp(get(c(i),'string'), 'Next shot'))
         nextb = c(i);
      elseif( strcmp(get(c(i),'string'), 'Previous shot'))
         prevb = c(i);
      elseif( strcmp(get(c(i),'string'), 'Done') )
         doneb = c(i);
      end
   end
end
if(isempty(rangeb))
   error('avgcvpinfo: did not find set range button!');
end
 
if(isempty(nextb))
   error('avgcvpinfo: did not find next pair button!');
end
 
if(isempty(prevb))
   error('avgcvpinfo: did not find previous pair button!');
end
if(isempty(doneb))
   error('avgcvpinfo: did not find done button!');
end
if( strcmp(action, 'get') )
   if( nargout < 5 )
      error('avgcvpinfo: not enough output arguements');
   else
      shotrange = get(doneb, 'userdata');
      shotlist = get(rangeb, 'userdata');
      tmp = get(nextb, 'userdata');
      if(~isempty(tmp))
         currentshot = tmp(1);
         stddev = tmp(2);
      else
         currentshot = [];
         stddev = [];
      end
      axeslist = get(prevb, 'userdata');
   end
end
if( strcmp(action, 'set') )
   if( nargin < 6 )
      error('avgcvpinfo: not enough input arguements');
   else
      set(rangeb, 'userdata', shotlist)
      if(isempty(currentshot))
         currentshot = 1;
      end
      if(isempty(stddev))
         stddev = 1;
      end
      set(nextb, 'userdata', [currentshot stddev]);
      set(prevb, 'userdata', axeslist);
      set(doneb, 'userdata', shotrange);
   end
end
