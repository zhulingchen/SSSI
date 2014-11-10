% TESTRAY: demo the v(x,z) raytrace code
%ray test
dg=10;

tstep=0:.004:3;
x=0:dg:5000;z=x';
v=1800+.6*(z*ones(1,length(x)));

rayvelmod(v,dg);

theta=pi*45/180;
r0=[0,0,sin(theta)/1800,cos(theta)/1800]';

tic
[t,r]=shootrayvxz(tstep,r0);
toc

tic
[tg,rg]=shootrayvxz_g(tstep,r0);
toc

tic
[tm,rm]=ode45('drayvec',tstep,r0);
toc

plot(rm(:,1),rm(:,2)+2*dg,rg(:,1),rg(:,2)+dg,r(:,1),r(:,2));flipy
legend('drayvec','shootrayvxz\_g','shootrayvxz')
