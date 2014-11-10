function dispavpstd(cvpavg, cvpstd)
f=gcf;
figure('menubar','none')
[nshots tmp] = size(cvpavg);
shots = 1:nshots;
% plot(cvpstd(:,1));
% hold on;
% plot(cvpavg(:,1));
errorbar(cvpavg(:,1), shots, cvpstd(:,1));
xlabel('shot number');
ylabel('xcoordinate');
title('Standard deviation with Cross over point average for the left hand side');
figure
plot(cvpstd(:,2));
hold on;
plot(cvpavg(:,2));
xlabel('shot number');
ylabel('xcoordinate');
title('Standard deviation with Cross over point average for the right hand side');
figure(f); set(gcf,'menubar','none');
