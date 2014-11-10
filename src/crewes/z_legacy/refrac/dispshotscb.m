function dispshotscb(action)
% Function determining the parameter for the arrival time curves display
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Display of which refracted arrivals (group per shot)?',...
             'Enter the first and last shot if "all" is not selected:',...
             'Enter the shot increment:');
   a=str2mat('all|shot (first,last)','10 25','5');
   askthingsinit('dispshotscb(''answer'')',q,a,[1 0 1],'Parameter for the refracted arrivals display');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   if( strcmp(deblank(a(1,:)), 'all') )
	shrange = 0
   else
	shrange = 1
   end
   sh = sscanf(a(2,:), '%d %d')
   sh1 = sh(1);
   sh2 = sh(2);
   inc = str2num(a(3,:));
   dispshots(shrange,sh1,sh2,inc);
end
