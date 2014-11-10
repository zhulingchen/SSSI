% RECSTUDY: Study record length limit on scattering angle
%
% Just run the script
ls{1}='-';ls{2}='-.';ls{3}=':';ls{4}='--';
figure;
vo=1800;c=.6;
T=[2. 4.  6. 8];
z=0:25:26000;
vo1=3500;c1=0;
for k=1:length(T)
	thetac=threc(T(k),vo1,c1,z);
    theta=threc(T(k),vo,c,z);
    ind=find(imag(theta)~=0);
    theta(ind)=nan*ind;
    ind=find(imag(thetac)~=0);
    thetac(ind)=nan*ind;
    h(2*k-1)=plot(z,thetac,['b' ls{k}]);
    if(k==1) hold; end
    h(2*k)=plot(z,theta,['r' ls{k}]);
end
set(gca,'ytick',[0:20:150])
ylabel('scattering angle in degrees')
xlabel('depth in meters')
legend(h,['Const v T=' int2str(T(1))],['Linear v(z) T=' int2str(T(1))],...
    ['Const v T=' int2str(T(2))],['Linear v(z) T=' int2str(T(2))],....
    ['Const v T=' int2str(T(3))],['Linear v(z) T=' int2str(T(3))],....
    ['Const v T=' int2str(T(4))],['Linear v(z) T=' int2str(T(4))])
%whitefig
%grid

	
