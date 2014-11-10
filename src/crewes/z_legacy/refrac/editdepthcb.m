function editdepthcb(action)
% Determine the parameter for the depth editing function
% Median filter window length
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Apply a median filter?',...
             'Enter the median filter window length');
   a=str2mat('Yes|No','5');
   askthingsinit('editdepthcb(''answer'')',q,a,[1 0 ],'Edit depth');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   if( strcmp(deblank(a(1,:)), 'Yes'))
	edit=0;
   else
	edit=1;
   end
   win=str2num(a(2,:));
   editdepth(edit,win);
end
