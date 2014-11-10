% Display of the arrival time curves sorted by shot
function dispshots(shrange,sh1,sh2,inc)
f=gcf;
fbtime=refdata('get','fbtime');
fbcoord=refdata('get','fbcoord');
nshots = refdata('get','nshots')
shotcoord = refdata('get','shotcoord');
% Arrival time curves for all incremental shot 
if (shrange==0)
   figcent(.6,.5)
   hold on;
   for n=1:inc:nshots
     plot(fbcoord(n,:),fbtime(n,:))
   end
   xy=axis;
   t=xy(4)-xy(3);
   d=t/40;
   for n=10:10:nshots         % Label every 10th shot
      str=sprintf('%d',n); 
      text(shotcoord(n),xy(3)+d,str)
   end
   text(xy(1)+100,xy(3)+2.5*d,'shot number')
   xlabel('Coordinate (m)');
   ylabel('Traveltime (ms)');
   title('Refracted arrivals');
   set(gcf,'units','pixels','menubar','none');
else
   % Arrival time curves for incremental shot in between two specified shot
   figcent(.6,.5);
   hold on;
   for n=sh1:inc:sh2
     plot(fbcoord(n,:),fbtime(n,:))
   end
   xy=axis;
   t=xy(4)-xy(3);
   d=t/10;
   for n=sh1:inc:sh2
     str=sprintf('%d',n); 
     text(shotcoord(n),xy(3)+d,str)
   end
   text(xy(1)+100,xy(3)+2*d,'shot number')
   xlabel('Coordinate (m)');
   ylabel('Traveltime (ms)');
   title('Refracted arrivals');
   set(gcf,'units','pixels','menubar','none');
end
figure(f);
