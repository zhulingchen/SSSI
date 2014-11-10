function dispcvpstat
% Display of the Cross Over point averages (left and rigth) with their 
% corresponding standard deviation and fold for each shot
% The Cross Over point averages are displayed in term of their coordinates (m) 
% and the interpolated time (ms) on the arrival curve 
f = gcf;
cvpi = refdata('get', 'cvpi');
cvpj = refdata('get', 'cvpj');
fbtime = refdata('get','fbtime');
fbcoord = refdata('get','fbcoord');
nshots = refdata('get', 'nshots');
shotcoord =refdata('get','shotcoord');
% Call of the Cross Over point averaging function
[cvpavg, cvpstd, cvpfold] = avgcvp(cvpi, cvpj, nshots);
cvpavg = refdata('get','cvpavg');
% Figure split in three part, first the Cross over point average values,
% second the corresponding standard deviation, and third the fold
figure('menubar','none');
set(gcf,'units','pixels','position',[0 0 864 720]);
coordaxes = axes('position',[.1 .7 .8 .25]);
stdaxes = axes('position', [.1 .4 .8 .25]);
foldaxes = axes('position',[.1 .1 .8 .25]);
% Cross Over point average values
axes(coordaxes);
hold on;
for n=1:nshots
   if (isnan(cvpavg(n,1)) ~=1)
	%timei=interp1(fbcoord(n,:),fbtime(n,:),cvpavg(n,1));
	%plot(cvpavg(n,1),timei,'color','c','linestyle','o')
	offset=shotcoord(n)-cvpavg(n,1);
	plot(shotcoord(n),offset,'color','c','linestyle','o')
   end
   if (isnan(cvpavg(n,2)) ~=1)
 	%timej=interp1(fbcoord(n,:),fbtime(n,:),cvpavg(n,2));
	%plot(cvpavg(n,2),timej,'color','r','linestyle','*')
	offset=cvpavg(n,2)-shotcoord(n);
	plot(shotcoord(n),offset,'color','r','linestyle','*')
   end
end
xy=axis;
t=xy(4)-xy(3);
d=t/10;
for n=10:10:nshots         % Label every 10th shot
    str=sprintf('%d',n); 
    text(shotcoord(n),xy(3)+d,str)
end
text(xy(1)+100,xy(3)+d,'shot number');
ylabel('Offset (m)')
title('Cross over point offset average from each shot for the left side (blue circle) and for the  right side (red star)')
% Standard deviation
axes(stdaxes);
hold on;
for n=1:nshots
    if (isnan(cvpavg(n,1)) ~=1)
	plot(shotcoord(n),cvpstd(n,1),'color','c','linestyle','o')
   end
   if (isnan(cvpavg(n,2)) ~=1)
	plot(shotcoord(n),cvpstd(n,2),'color','r','linestyle','*')
   end
end
xy=axis;
t=xy(4)-xy(3);
d=t/10;
for n=10:10:nshots         % Label every 10th shot
    str=sprintf('%d',n); 
    text(shotcoord(n),xy(3)+d,str)
end
text(xy(1)+100,xy(3)+d,'shot number');
ylabel('Standard deviation (m)')
% Fold
axes(foldaxes);
hold on;
for n=1:nshots
    if (isnan(cvpavg(n,1)) ~=1)
	plot(shotcoord(n),cvpfold(n,1),'color','c','linestyle','o')
   end
   if (isnan(cvpavg(n,2)) ~=1)
	plot(shotcoord(n),cvpfold(n,2),'color','r','linestyle','*')
   end
end
xy=axis;
t=xy(4)-xy(3);
d=t/10;
for n=10:10:nshots         % Label every 10th shot
    str=sprintf('%d',n); 
    text(shotcoord(n),xy(3)+d,str)
end
text(xy(1)+100,xy(3)+d,'shot number');
ylabel('Fold')
xlabel('Coordinate (m)')
%set(gcf,'units','pixels','position',[0 0 864 720]);
figure(f); set(gcf,'menubar','none');
