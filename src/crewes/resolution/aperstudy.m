% APERSTUDY: study the aperture effect
%
% Just run the script
ls{1}='-';ls{2}='-.';ls{3}=':';ls{4}='--';
figure;
vo=1800;c=.6;
A=[1000  4000  12000  20000];
z=0:25:20000;
vo1=3500;c1=0;
h=zeros(1,8);
for k=1:length(A)
    theta=thaper(A(k),vo,c,z);
	thetac=thaper(A(k),vo1,c1,z);
    h(2*k-1)=plot(z,thetac,['b' ls{k}]);
    if(k==1) hold; end
    h(2*k)=plot(z,theta,['r' ls{k}]);
end




%p=get(gcf,'position');
%set(gcf,'position',[p(1:2) 700 700])
ylabel('scattering angle in degrees')
xlabel('depth in meters')
set(gca,'xtick',[0:5000:20000])
set(gca,'ytick',[0:30:180])
set(gca,'ylim',[0 180])
legend(h,['Const v A=' int2str(A(1))],['Linear v(z) A=' int2str(A(1))],...
    ['Const v A=' int2str(A(2))],['Linear v(z) A=' int2str(A(2))],....
    ['Const v A=' int2str(A(3))],['Linear v(z) A=' int2str(A(3))],....
    ['Const v A=' int2str(A(4))],['Linear v(z) A=' int2str(A(4))])
%whitefig
%grid
	
