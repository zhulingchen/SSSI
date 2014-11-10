function [hg1]=depthreject(hg1,dev,standard)
[coord, alldepth, depth, fold] = avgdepth(fbcoord,hg1);
for n=1:nshots
	% Average the crossover points for 'i' (left side)
	validcvpi = find(~isnan(cvpi(n,:)));
	x=cvpi(n,validcvpi);
	d = abs(x - meani(n));
    if (standard==1)
	f = dev * devi(n);
    else
	f = dev;
    end
        badcvpi = find(d>f);
	[a b]=size(badcvpi);
	if (b ~=0)
	   x(badcvpi)=NaN*badcvpi;
	   cvpi(n,validcvpi) = x;
	end
	% Average the crossover points for 'j' (right side)
	validcvpj = find(~isnan(cvpj(n,:)));
	x=cvpj(n,validcvpj);
	d = abs(x - meanj(n));
    if (standard==1)
	f = dev * devj(n);
    else
	f = dev;
    end
        badcvpj = find(d>f);
	[a b]=size(badcvpj);
	if (b ~=0)
	   x(badcvpj)=NaN*badcvpj;
	   cvpj(n,validcvpj) = x;
	end
end
