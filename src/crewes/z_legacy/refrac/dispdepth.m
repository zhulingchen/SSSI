function dispdepth
% Display of the calculated depth model and the corresponding 
% standard deviation and fold for each receiver  
f=gcf;
depth=refdata('get','depth');
recelev=refdata('get','recelev');
% Find the position of the layer1-layer2 interface according 
% to the surface elevation
depthelev=recelev(2,:)-depth(2,:);
% Figure split in three part, first the surface elevation and depth average values,
% second the corresponding standard deviation, and third the fold
figure('menubar','none');
coordaxes = axes('position',[.1 .7 .8 .25]);
stdaxes = axes('position', [.1 .4 .8 .25]);
foldaxes = axes('position',[.1 .1 .8 .25]);
% Depth average values with surface elevation
axes(coordaxes);
hold on;
plot(depth(1,:),recelev(2,:),'color','y','linestyle','-')
plot(depth(1,:),depthelev,'color','g','linestyle','-.')
xmin=min(depth(1,:));
xmax=max(depth(1,:));
ymin=min(depthelev(1,:))-20;
ymax=max(recelev(2,:))+20;
axis([xmin xmax ymin ymax]);
ylabel('Elevation (m)')
title('Depth model (surface= -yellow; layer 1-2 interface= -.green)')
% Standard deviation
axes(stdaxes);
hold on;
plot(depth(1,:),depth(4,:),'color','g','linestyle','*')
ylabel('Standard deviation (m)')
% Fold
axes(foldaxes);
hold on;
plot(depth(1,:),depth(3,:),'color','g','linestyle','+')
ylabel('Fold')
xlabel('Coordinate (m)')
set(gcf,'units','pixels','position',[0 0 864 720]);
figure(f); set(gcf,'menubar','none');
