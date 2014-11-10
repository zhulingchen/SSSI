function rangecb(action)
if( nargin < 1 )
   action = 'init';
end
[shots, k, axeslist] = cvpinfo('get');
if( strcmp(action,'init'))
   q=str2mat('Enter the shot pair range: (first left:last left),(first right:last right)');
   a=str2mat('(10:15),(11:20)');
   askthingsinit('rangecb(''answ'')',q,a,[1],'Shot pair range');
elseif( strcmp(action,'answ'))
   a=askthingsfini;
   [strings tmp] = size(a);
   answer=a(1,:)
   % Calculate an array of shot pairs (use the answer from the popup box)
   shots = shotrange(answer);
   k = 1;
   % Turn on the next and previous buttons
   c = get(gcf, 'children');
   [nc tmp] = size(c);
   for i=1:nc
      if( strcmp(get(c(i),'type'),'uicontrol') )
         if( strcmp( get(c(i), 'string'), 'Next pair' ) | ... 
             strcmp( get(c(i), 'string'), 'Previous pair' ) )
            set(c(i), 'visible', 'on');
         end
      end
   end
   cvpinfo('set',shots, k, axeslist);
   editcvp('run')
end
