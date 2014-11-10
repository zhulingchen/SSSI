function avgnextcb(action)
[shotrange shotINClist currentshot axeslist stddev] = avgcvpinfo('get');
nshots = refdata('get','nshots');
if( strcmp(action,'next'))
   if(isempty(currentshot))
      currentshot = 1;
   else
      currentshot = currentshot + 1;
   end
end
if( strcmp(action,'prev'))
   if(isempty(currentshot))
      currentshot = nshots;
   else
      currentshot = currentshot - 1;
   end
end
if( currentshot < 1 ) 
   currentshot = 1;
   disp('First shot reached.');
end
 
%numshotpairs = length(shotrange);
if( currentshot >  nshots )
   currentshot = nshots;
   disp('End of shots reached.');
end
avgcvpinfo('set', shotrange, shotINClist, currentshot, axeslist, stddev);
editcvpavg('run');
