function disppoly(action)
% Determination of the parameter for the polyfit of the
% cross over point averages
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Enter the polynomial degree:');
 a=str2mat('5');
   askthingsinit('disppoly(''answer'')',q,a,[1],'Parameter for the CVP average Polyfit function');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   cvpi=refdata('get','cvpi');
   cvpj=refdata('get','cvpj');
   nshots=refdata('get','nshots');
   shotcoord=refdata('get','shotcoord');
   % Call the CVP averaging function
   [cvpavg cvpstd cvpfold]=avgcvp(cvpi,cvpj,nshots);
   f=gcf;
      % Plot the CVP averages with the Polyfit curves
      poly=str2num(a(1,:));
      goodleft=find(~isnan(cvpavg(:,1)));
      goodright=find(~isnan(cvpavg(:,2)));
      offsetleft=shotcoord(goodleft)-cvpavg(goodleft,1);
      offsetright=cvpavg(goodright,2)-shotcoord(goodright);
      l=polyfit(goodleft,offsetleft,poly);
      r=polyfit(goodright,offsetright,poly);
      left=polyval(l,goodleft);
      right=polyval(r,goodright);
      figure('menubar','none')
      hold on
      plot(goodleft,offsetleft,'co')
      plot(goodleft,left,'c-.')
      plot(goodright,offsetright,'r*')     
      plot(goodright,right,'r')
      title('Cross over point offset average from each shot for the left side (blue circle) and for the right side (red star)')
      ylabel('Offset(m)')
      xlabel('Shot number')
      set(gcf,'units','pixels','position',[0 0 864 576]);
      figure(f)
end
