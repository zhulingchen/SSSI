load demodata
plotimage(seisfilt,t,x);
title('Time Section (Explodogram)')
whitefig;bigfont(gca,1.5,1);
ylabel('seconds');xlabel('meters')
set(gcf,'position',[463   306   560   420])
plotimage(r,z,x);
title('Depth Model');
whitefig;bigfont(gca,1.5,1);
ylabel('meters');xlabel('meters')
set(gcf,'position',[2   307   560   420])
rayvelmod(vel,dx)
figure(1)
