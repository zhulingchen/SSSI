load demodata2
plotimage(seisfilt,t,x);
title('Time Section (Explodogram)')
whitefig;bigfont(gca,1.5,1);
ylabel('seconds')
set(gcf,'position',[463   306   560   420])
plotimage(vel-.5*(vhigh+vlow),z,x)
title('Depth Model');
whitefig;bigfont(gca,1.5,1);
ylabel('depth')
set(gcf,'position',[2   307   560   420])
rayvelmod(vel,dx)
figure(1)
