function editvelcb(action)
% Determine the parameter for the velocity editing function
% median filter window length
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Apply a median filter ?',...
	'To which velocities?',...
	'Median filter window length (odd number):');
   
   a=str2mat('Yes|No','First|Second|Both','5');
   askthingsinit('editvelcb(''answer'')',q,a,[1 0 0],'Edit velocity');
   
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   if( a ~= -1 )
     [strings tmp] = size(a);
     if( strcmp(deblank(a(1,:)), 'Yes'))
	edit=0;
     else
	edit=1;
     end
     
     if( strcmp(deblank(a(2,:)), 'First'))
	v=1;
     elseif( strcmp(deblank(a(2,:)), 'Second'))
	v=2;
     else
	v=3;
     end
     
     win=str2num(a(3,:));
     
     editvel(edit,v,win);
   end
end
