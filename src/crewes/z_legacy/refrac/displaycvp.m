figure('menubar','none')
hold on;
plot(fbcoord(j,:),fbtime(j,:))
plot(fbcoord(i,:),fbtime(i,:))
cvpi(i,j)
cvpj(i,j)
if (~isnan(cvpi(i,j)))
	indexi=find(fbcoord(i,:)==cvpi(i,j));
	timei=fbtime(i,indexi);
	plot(cvpi(i,j),timei,'o')
end
if (~isnan(cvpj(i,j)))
	indexj=find(fbcoord(j,:)==cvpj(i,j));
	timej=fbtime(j,indexj);
	plot(cvpj(i,j),timej,'o')
end
