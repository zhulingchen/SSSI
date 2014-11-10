function dispstatic
% Display receiver and shot static corrections for the weathering layer
% (first layer, the elevation and the total
f=gcf;
recstat=refdata('get','recstat');
shotstat=refdata('get','shotstat');
recelev=refdata('get','recelev');
shotcoord=refdata('get','shotcoord');
nshots=refdata('get','nshots');
nrecs=refdata('get','nrecs');
% Receiver static corrections
figure('menubar','none');
hold on;
plot(recstat(1,:),'color','c','linestyle','+')
plot(recstat(2,:),'color','r','linestyle','*')
plot(recstat(3,:),'color','g','linestyle','o')
xlabel('Receiver number')
ylabel('Static corrections (ms)')
title('Receiver static corrections (weathering=plus; elevation=star; total=circle)')
% Shot static corrections
figure('menubar','none');
hold on;
plot(shotstat(1,:),'color','c','linestyle','+')
plot(shotstat(2,:),'color','r','linestyle','*')
plot(shotstat(3,:),'color','g','linestyle','o')
xlabel('Shot number')
ylabel('Static corrections (ms)')
title('Shot static corrections (weathering=plus; elevation=star; total=circle)')
% Surface consistent statics (X-coordinate)
figure('menubar','none');
hold on;
plot(recelev(1,:),recstat(1,:),'color','c','linestyle','+')
plot(recelev(1,:),recstat(2,:),'color','r','linestyle','*')
plot(recelev(1,:),recstat(3,:),'color','g','linestyle','o')
xlabel('Coordinate (m)')
ylabel('Static corrections (ms)')
title('Surface consistent static corrections (weathering=plus; elevation=star; total=circle)')
xy=axis;
t=xy(4)-xy(3);
d=t/20;
for n=10:10:nshots         % Label every 10th shot
    str=sprintf('%d',n); 
    text(shotcoord(n),xy(3)+d,str)
end
text(xy(1)+100,xy(3)+2*d,'shot number')
for n=20:20:nrecs         % Label every 20th receiver
    str=sprintf('%d',n); 
    text(recelev(1,n),xy(3)+3*d,str)
end
text(xy(1)+100,xy(3)+4*d,'receiver number')
set(gcf,'units','pixels','position',[0 0 864 720],'menubar','none');
figure(f);
