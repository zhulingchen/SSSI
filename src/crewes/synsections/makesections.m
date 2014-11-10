% MAKESECTIONS: demo the use of the section tools
%
% Just run the script

dt=.004;dx=10;tmax=2;xmax=2500;v=3000;
[w,tw]=ricker(.004,40);
tic
[amat,t,x]=makestdsyn(dt,dx,tmax,xmax,v,w,tw);
toc
tic
[amath,t,x]=makestdsynh(dt,dx,tmax,xmax,v,w,tw);
toc
plotimage(amat,t,x);
title('Made without diffractions')
[fk,f,k]=fktran(amat,t,x);
plotimage(abs(fk),f,k);
xlabel('wavenumber');ylabel('Hertz')
title('FK spectrum without diffractions')

plotimage(amath,t,x);
title('Made with diffractions')
[fkh,f,k]=fktran(amath,t,x);
plotimage(abs(fkh),f,k);
xlabel('wavenumber');ylabel('Hertz')
title('FK spectrum with diffractions')

%velocity model for raytracing
zmax=v*max(t)/2;
z=0:dx:zmax;

vmodel=v*ones(length(z),length(x)); %build model
rayvelmod(vmodel,dx); % initialize for raytracing
