function avgdepthcb(action)
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Use a constant velocity for the first layer?',...
             'Enter the constant velocity (m/s):');
   a=str2mat('No|Yes','650')
   askthingsinit('avgdepthcb(''answer'')',q,a,[1 0],'Parameter for the depth calculation');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   plust=refdata('get','plust');
   fbcoord =refdata('get','fbcoord');
   td1 = refdata('get', 'td1');
   v2rec = refdata('get','v2rec');
   if(strcmp( deblank(a(1,:)),'No'))
     v1rec = refdata('get','v1rec');
   else
     v1 = str2num(a(2,:));
     v1 = v1/1000;
     [m n] = size(v2rec);
     v1rec = v1*ones(1,n);
   end
   [depth] = calcdepth(plust,v1rec,v2rec);
   refdata('set','depth',depth);
   % Update menus
   PMTsetmenus;
end
