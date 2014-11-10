function staticcb(action)
% Determine the parameter for the static function
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Enter datum elevation (m):',...
             'Use a pseudo-datum?',...
             'Use the 2nd layer velocity as replacement velocity?',...
             'Replacement velocity (m/s):',...
             'Use a constant velocity for the first layer?',...
             'Constant velocity (m/s):');
   a=str2mat('100','Yes|No','No|Yes','1600','No|Yes','600');
   askthingsinit('staticcb(''answ'')',q,a,[1 1 1 0 1 0],'Static computation');
elseif( strcmp(action,'answ'))
   a=askthingsfini;
   [strings tmp] = size(a);
   datum=str2num(a(1,:));
   repvel=str2num(a(4,:));
   if( strcmp( deblank(a(2,:)), 'Yes') )
	psd=1; % use of a pseudo-datum (method 2)
   else
	psd=0; % do not use a pseudo-datum (method 1)
   end
   if( strcmp( deblank(a(3,:)), 'Yes') )
	v2rep=1; % use of a v2 as a replacement velocity
   else
	v2rep=0; % use of a constant for replacement velocity
   end
   if(strcmp( deblank(a(5,:)),'No'))
     % Use of the first layer velocity as calculated
     v1rec = refdata('get','v1rec'); 
   else
     % Use of a specified constant first layer velocity
     v1 = str2num(a(6,:));
     [m n] = size(v2rec);
     v1rec = v1*ones(1,n);
   end
   v2rec = refdata('get','v2rec');
   depth=refdata('get','depth');
   shotcoord=refdata('get','shotcoord');
   shotelev=refdata('get','shotelev');
   recelev=refdata('get','recelev');
   uphole=refdata('get','uphole');
   % Call the static function
   [recstat,shotstat] =static(depth, shotcoord, shotelev, recelev, uphole, datum,repvel,psd, v1rec, v2rec, v2rep);
   refdata('set','recstat',recstat);
   refdata('set','shotstat',shotstat);
end
