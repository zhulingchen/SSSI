function disprectime(rtrange,rt1,mint);
% Display of the reciprocal time for shot pair(s) 
f=gcf;
diffmat=refdata('get','diffmat');
nshots=refdata('get','nshots');
absdiffmat=abs(diffmat);
% All the possible shot pairs
if (rtrange==0)
	figure('menubar','none')
	image(absdiffmat);
	set(gca,'ydir','normal');
	colormap(hsv);
        colorbar;
	xlabel('shot number');
	ylabel('shot number');
	zlabel('Reciprocal traveltime difference (ms)');
	title('Reciprocal times ')
	figure('menubar','none');
	[s1,s2]=find(absdiffmat>mint);
	plot(s2,s1,'*');
	xlabel('shot number');
	ylabel('shot number');
	titlestr = sprintf('Shot pairs with over %d ms of reciprocal time difference ', mint);
	title(titlestr);
% Only two shots with all the other possible shots 
else
        figure('menubar','none')
	hold on;
	plot(1:nshots,absdiffmat(rt1,:),'gx')
	plot(1:nshots,absdiffmat(:,rt1),'gx')
	a1=find(absdiffmat(rt1,:)>mint);
	a2=find(absdiffmat(:,rt1)>mint);
	plot(a1,absdiffmat(rt1,a1),'go');
	plot(a2,absdiffmat(a2,rt1),'go');
	xlabel('shot number');
	ylabel('Reciprocal traveltime difference (ms)');
	titlestr = sprintf('Reciprocal time difference for shot %d (circles show a time difference over %d) ',rt1, mint);
	title(titlestr);
        set(gcf,'units','pixels','position',[0 0 864 576])
end
