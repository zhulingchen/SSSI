function avgrangecb(action)
 
if( nargin < 1 )
   action = 'init';
end
 
[shotrange shotINClist shot axeslist stddev] = avgcvpinfo('get');
 
if( strcmp(action,'init'))
   nshots = refdata('get', 'nshots');
   q = str2mat('First shot number:', 'Last shot number:');
   if(~isempty(shotrange))
      rangeend = max(shotrange);
      if( rangeend > nshots )
         rangeend = nshots;
      end
      a = str2mat( num2str(shotrange(1)), num2str(rangeend) );
   else
      a = str2mat( '1', sprintf('%d',nshots));
   end
   askthingsinit('avgrangecb(''answ'')',q,a,[1 1],'Shot range');
elseif( strcmp(action, 'answ') )
   a = askthingsfini;
   % a is -1 if 'cancel' was pressed
   if( a ~= -1 )
      [strings tmp] = size(a);
 
      shotrange = str2num(a(1,:)) : str2num(a(2,:));
      shot = shotrange(1);
   
      if(isempty(shotINClist))
         nshots = refdata('get', 'nshots');
         shotINClist = 1:nshots;
      end
      avgcvpinfo('set', shotrange, shotINClist, shot, axeslist, stddev);
      % Turn on the buttons that can now work.
      c = get(gcf, 'children');
      [nc tmp] = size(c);
      for i=1:nc
         if( strcmp(get(c(i),'type'),'uicontrol') )
            if( strcmp(get(c(i),'string'), 'Next shot') | ...
               strcmp(get(c(i),'string'), 'Previous shot') | ...
               strcmp(get(c(i),'string'), 'Recompute') )
               set(c(i), 'visible', 'on');
            end
         end
      end
      editcvpavg('label');
   end
end
  
