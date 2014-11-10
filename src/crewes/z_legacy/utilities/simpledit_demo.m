% simpledit_demo
%generate two spiky reflectivities
[r1,t]=reflec(1.0,.004,.5,3,pi);
[r2,t]=reflec(1.0,.004,.5,3,exp(1));
%make two fake horizons
x=1:.5:10;
y1=.4+(x-1)*(-.033);
y2= .8 -.4*cos( (x-4.5)/9 );
%plot
figure;
p=get(gcf,'position');
set(gcf,'position',[p(1:2) 600 p(4)]);
h1=line( x,y1, 'color','r');
h2=line( x,y2, 'color','c');
h3=line( r1+3, t, 'color','b');
h4=line( r2+6, t, 'color','b');
set(gca,'ydir','reverse');
simpledit
editlines('multnanoff');
