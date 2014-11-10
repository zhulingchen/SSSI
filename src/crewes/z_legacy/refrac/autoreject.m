function [cvpi,cvpj]=autoreject(cvpi,cvpj,nshots,dev,standard)
% Function used to reject some of the cross over point based on the standard
% deviation or a constant difference from the cross over point average for
% each shot
si=NaN.*ones(1,nshots);
sj=NaN.*ones(1,nshots);
mi=NaN.*ones(1,nshots);
mj=NaN.*ones(1,nshots);
% Call the function calculating the average cross over points 
[cvpavg cvpstd,cvpfold] = avgcvp(cvpi, cvpj,nshots);
% Cross over point average for each shot (cvpi=left side; cvpj=right side)
meani = cvpavg(:,1);
meanj = cvpavg(:,2);
devi = cvpstd(:,1);
devj = cvpstd(:,2);
% Determination of the bad cross over point based on the standard deviation 
% or a constant limit
for n=1:nshots
	% Average the cross over points for 'cvpi' (left side)
	validcvpi = find(~isnan(cvpi(n,:)));
	x=cvpi(n,validcvpi);
	d = abs(x - meani(n));
    if (standard==1)
	f = dev * devi(n); % rejection based on the standard deviation
    else
	f = dev; % rejection based on a constant limit
    end
        badcvpi = find(d>f);
	[a b]=size(badcvpi);
	if (b ~=0)
	   x(badcvpi)=NaN*badcvpi;
	   cvpi(n,validcvpi) = x;
	end
	% Average the cross over points for 'cvpj' (right side)
	validcvpj = find(~isnan(cvpj(:,n)));
	x=cvpj(validcvpj,n);
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
	   cvpj(validcvpj,n) = x;
	end
end
