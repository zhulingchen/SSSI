function nextcb(action)
[shots, k, axeslist] = cvpinfo('get');
[numshotpairs tmp] = size(shots);
if( strcmp(action,'next') )
   if(isempty(k))
      k = 1;
   else
      k = k+1;
   end
end
if( strcmp(action,'prev') )
   if(isempty(k))
      k = numshotpairs;
   else
      k = k-1;
   end
end
if( k < 1 ) 
   k = 1;
   disp('Beginning of shot pair range reached.');
end
if( k > numshotpairs )
   k = numshotpairs;
   disp('End of shot pair range reached.');
end
cvpinfo('set', shots, k, axeslist);
editcvp('run');
