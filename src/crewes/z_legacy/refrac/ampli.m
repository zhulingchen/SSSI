absseisp=abs(seisp);
absseisr=abs(seisr);
trcz=near(seisp(1,:),trc1,trc2);
timez=near(seisp(:,1),t1,t2);
figure
plot(absseisr(timez,trcz),absseisp(timez,trcz),'o');
xlabel('Radial Component Amplitude')
ylabel('Vertical Component Amplitude')
titlestr = sprintf('Relative Amplitude for source gather #1 trace number %d to %d between %d and %d ms',trc1,trc2,t1,t2);
title(titlestr);
axis('equal')
axis('square')
