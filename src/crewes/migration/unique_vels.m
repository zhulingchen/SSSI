function vrefs=unique_vels(model)
[rm,cm]=size(model);
N=max(rm,cm);
vs=zeros(1,N);
count=0;
for j=1:N;
	temp=model;
	if sum(temp)>0
		inds=find(model~=temp(1));
		model=temp(inds);
		vs(j)=temp(1);
		count=count+1;
	end
end
vrefs=vs(1:count);
