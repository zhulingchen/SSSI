figure('menubar','none')
hold on;
for n=1:101
	plot(fbcoord(n,:),fbtime(n,:))
    if (isnan(cvpavg(n,1)) ~=1)
	%indexi=find(fbcoord(n,:)==cvpavg(n,1));
	%timei=fbtime(n,indexi);
	timei=interp1(fbcoord(n,:),fbtime(n,:),cvpavg(n,1));
	plot(cvpavg(n,1),timei,'o')
    end
   if (isnan(cvpavg(n,2)) ~=1)
	%indexj=find(fbcoord(n,:)==cvpavg(n,2));
	%timej=fbtime(n,indexj);
 	timej=interp1(fbcoord(n,:),fbtime(n,:),cvpavg(n,2));
	plot(cvpavg(n,2),timej,'o')
    end
end
