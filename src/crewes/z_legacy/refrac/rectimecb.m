function rectimecb(action)
% Determine the parameter for the reciprocal time difference function
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Reciprocal time check on which shot pairs ?',...
             'Enter the shot pair if "all" is not selected:',...
             'Identify the shot pairs having a reciprocal time check difference over (ms)?');
   a=str2mat('all|shot pair (i,j)','10 15','5');
   askthingsinit('rectimecb(''answer'')',q,a,[1 0 1],'Parameter for the Reciprocal time check');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   if( strcmp(deblank(a(1,:)), 'all') )
      rtrange = 0;
   else
      rtrange = 1;
   end
   refdata('set','rtrange',rtrange);
   disprectimem = refdata('get','disprectimem');
   set(disprectimem,'enable','on');
   rtpair = sscanf(a(2,:), '%d %d')
   rtpair1 = rtpair(1);
   rtpair2 = rtpair(2);
   refdata('set', 'rtpair1', rtpair1 );
   refdata('set', 'rtpair2', rtpair2 );
   mint = str2num(a(3,:));
   refdata('set', 'mint', mint);
   fbcoord = refdata('get','fbcoord');
   shotcoord = refdata('get','shotcoord');
   fbtime = refdata('get','fbtime');
   nshots = refdata('get','nshots');
   % Call the reciprocal time difference function
   [diffmat]=rectime(rtrange,rtpair1,rtpair2,fbcoord,shotcoord,fbtime,mint,nshots);
 
   % Only save the diffmatrix if all shot pairs have been computed
   if( rtrange == 0 )
      refdata('set','diffmat',diffmat);
   end
end
