function cvpavgSRcb(action)
if( nargin < 1 ) 
   action = 'init';
end
nshots = refdata('get', 'nshots');
[shotrange shotlist shot axeslist stddev] = avgcvpinfo('get');
if( strcmp(action, 'init') )
   window = refdata('get','window');
   windowmn = refdata('get','windmn');
   nd = refdata('get', 'nd');
   sgap = refdata('get', 'shotgap');
   slen = refdata('get', 'shotlength');
   flatslope = refdata('get', 'flatslope');
   q = str2mat( 'Maximum "flat" slope (%):', ...
                'Standard deviation limit:', ...
                'Gap from shot:', ...
                'Number of  shots/side:', ...
                'Differentiation seperation:', ...
		'Median filter window:', ...
		'Mean filter window:' );
   a = str2mat( sprintf('%3.1f',flatslope*100), num2str(stddev), ...
                num2str(sgap), num2str(slen), num2str(nd), ...
                num2str(window), num2str(windowmn) );
   askthingsinit('cvpavgSRcb(''answer'')', q, a, [1 1 1 1 1 1 1], 'Parameters');
elseif( strcmp(action, 'answer') )
   a = askthingsfini;
   if( a ~= -1 )
      [strings tmp] = size(a);
      flatslope = str2num(deblank(a(1,:))) / 100;
      stddev = str2num(deblank(a(2,:)));
      sgap = str2num(deblank(a(3,:)));
      slen = str2num(deblank(a(4,:)));
      nd = str2num(deblank(a(5,:)));
      window = str2num(deblank(a(6,:)));
      windmn = str2num(deblank(a(7,:)));
      if(isempty(shotlist))
         shotlist = 1:nshots;
         shot = 1;
      else
         shot1 = max(1,shot-slen-sgap);
         shot2 = max(1,shot-sgap);
         shot3 = min(nshots, shot+sgap);
         shot4 = min(nshots, shot+sgap+slen);
         shotlist = [shot1:shot2 shot3:shot4];
      end
      avgcvpinfo('set', shotrange, shotlist, shot, axeslist, stddev);
      refdata('set', 'shotgap', sgap);
      refdata('set', 'shotlength', slen);
      refdata('set', 'flatslope', flatslope );
      refdata('set', 'nd', nd);
      refdata('set', 'window', window);
      refdata('set', 'windmn', windmn);
      editcvpavg('label');
   end
end
