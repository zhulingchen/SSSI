% ALIASTUDY: study aliasing versus depth
%
% Just run the script
ls{1}='-';ls{2}='-.';ls{3}=':';ls{4}='--';
figure;
vo=1800;c=.6;
dx=[20 40 80 160];
z=0:25:20000;
f=60;
vo1=3500;c=01;

for k=1:length(dx)
	thetac=thalias(dx(k),f,vo1,c1,z);
    theta=thalias(dx(k),f,vo,c,z);
    h(2*k-1)=plot(z,thetac,['b' ls{k}]);
    if(k==1) hold; end
    h(2*k)=plot(z,theta,['r' ls{k}]);
end

%p=get(gcf,'position');
%set(gcf,'position',[p(1:2) 700 700])
ylabel('scattering angle in degrees')
xlabel('depth in meters')
set(gca,'xtick',[0:2500:20000])
set(gca,'ytick',[0:30:90])
legend(h,['Const v \Delta x=' int2str(dx(1))],['Linear v(z) \Delta x=' int2str(dx(1))],...
    ['Const v \Delta x=' int2str(dx(2))],['Linear v(z) \Delta x=' int2str(dx(2))],....
    ['Const v \Delta x=' int2str(dx(3))],['Linear v(z) \Delta x=' int2str(dx(3))],....
    ['Const v \Delta x=' int2str(dx(4))],['Linear v(z) \Delta x=' int2str(dx(4))])
%whitefig
%grid
	
