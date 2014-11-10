function timerejectcb(action)
if( nargin < 1)
   action = 'init';
end
if( strcmp(action,'init'))
   q = str2mat('Autorejection based on constant or standard devation limit?',...
               'Enter the limit:');
   a = str2mat('Standard deviation|Constant limit','1');
   askthingsinit('timerejectcb(''answer'')',q,a,[1 1],'Autoreject time parameters');
elseif( strcmp(action, 'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   dev = str2num(a(2,:));
   refdata('set', 'dev', dev);
   if(strcmp( deblank(a(1,:)),'standard deviation'))
	standard=1;
   else
	standard=0;
   end
   refdata('set','standard',standard);
   fbcoord = refdata('get', 'fbcoord');
   td1 = refdata('get', 'td1');
   plustreject = 1;     % True, do the rejection
   plust = avgplustime(fbcoord, td1, standard, plustreject, dev);
   refdata('set','plust', plust);
end
