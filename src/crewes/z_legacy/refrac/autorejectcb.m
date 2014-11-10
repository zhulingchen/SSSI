function autorejectcb(action)
% Determination of the parameter for the cross over point rejection
% The rejection is based on a factor of the standard deviation or on a
% constant limit difference from the average cross over point for each shot  
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Constant or standard deviation limit?','Enter the limit:');
   a=str2mat('standard deviation|constant limit','1');
   askthingsinit('autorejectcb(''answer'')',q,a,[1 1],'Parameter for the Autorejection of CVP');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   dev = str2num(a(2,:));
   refdata('set', 'dev', dev);
   if(strcmp( deblank(a(1,:)),'standard deviation'))
	standard=1 % rejection based on the standard deviation
   else
	standard=0 % rejection based on a constant limit
   end
   refdata('set','standard',standard);
   cvpi = refdata('get', 'cvpi');
   cvpj = refdata('get', 'cvpj');
   nshots = refdata('get', 'nshots');
   [cvpi, cvpj] = autoreject(cvpi,cvpj,nshots,dev,standard);
   refdata('set','cvpi',cvpi);
   refdata('set','cvpj',cvpj);
end
