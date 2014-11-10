function [cvpavg,cvpstd,cvpfold]=avgcvp(cvpi,cvpj,nshots)
% Calculation of the average cross over points for each shot
% The left side cross over point (cvpi) and the rigth side cross over point
% (cvpj) are averaged independently and store in two collom
% The standard deviation and the fold are found for each of them
 
cvpavg=NaN.*ones(nshots,2);
cvpstd=NaN.*ones(nshots,2);
cvpfold=NaN.*ones(nshots,2);
% Averaging loop for each shot
for n=1:nshots
    % Left cross over points averaging
	validcvpi = find(~isnan(cvpi(n,:)));
	[a b]= size(validcvpi);
    if( b ~=0 )
	avgcvpi(n)=mean(cvpi(n,validcvpi));
	cvpavg(n,1)=avgcvpi(n);
        cvpstd(n,1)=std(cvpi(n,validcvpi));
	cvpfold(n,1)=length(validcvpi);
    end
    % Rigth cross over points averaging
	validcvpj = find(~isnan(cvpj(:,n)));
	[a b]= size(validcvpj);
    if( b ~=0 )
	avgcvpj(n)=mean(cvpj(validcvpj,n));
	cvpavg(n,2)=avgcvpj(n);
        cvpstd(n,2)=std(cvpj(validcvpj,n));
	cvpfold(n,2)=length(validcvpj);
    end
end
