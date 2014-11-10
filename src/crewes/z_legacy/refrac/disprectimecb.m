function disprectimecb(action)
% Function determining the parameter for the reciprocal time display
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Display which reciprocal times (shot group)?',...
             'Enter the shot group if "all" is not selected:',...
             'Identify shot pairs over a minimum time difference?');
   a=str2mat('all|shot (i,j)','25','5');
   askthingsinit('disprectimecb(''answer'')',q,a,[1 0 1],'Parameter for the reciprocal times display');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   if( strcmp(deblank(a(1,:)), 'all') )
	rtrange = 0
   else
	rtrange = 1
   end
   rt1 = str2num(a(2,:));
   mint = str2num(a(3,:));
   disprectime(rtrange,rt1,mint);
end
