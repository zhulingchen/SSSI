function dispplust
% Display of the average Plus Time values with their
% corresponding standard deviation and fold for each receiver
f=gcf;
fbcoord=refdata('get','fbcoord');
plust=refdata('get','plust');
% Figure split in three part, first the Plus Time average values,
% second the corresponding standard deviation, and third the fold
figure('menubar','none');
coordaxes = axes('position',[.1 .7 .8 .25]);
stdaxes = axes('position', [.1 .4 .8 .25]);
foldaxes = axes('position',[.1 .1 .8 .25]);
% Plus Time average values
axes(coordaxes);
hold on;
plot(plust(1,:),plust(2,:),'color','b','linestyle','-.')
ylabel('Plus Time (ms)')
title('Plus time')
% Standard deviation
axes(stdaxes);
hold on;
plot(plust(1,:),plust(4,:),'color','b','linestyle','*')
ylabel('Standard deviation (ms)')
% Fold
axes(foldaxes);
hold on;
plot(plust(1,:),plust(3,:),'color','b','linestyle','+')
ylabel('Fold')
xlabel('Coordinate (m)')
set(gcf,'units','pixels','position',[0 0 864 720]);
figure(f); set(gcf,'menubar','none');
