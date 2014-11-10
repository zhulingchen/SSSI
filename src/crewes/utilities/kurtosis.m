function k=kurtosis(X,flag)
% k=kurtosis(X), k=kurtosis(X,flag)
%
% Kurtosis is a measure of the peakedness of the probablility distribution
% of a real random variable.  It is used in proabablitiy and statistics.
%
% This function will return the kurtosis of set of data.  If the data is a
% 2-D matrix it will return a value for each column.
%
% if flag=0 a scaling is applied. flag=1; is the default.
%
if nargin <2
    flag=1;
end


sz=size(X);
if sz(1)==1
    X=X';
sz=size(X);
end
n=sz(1);

xbar=mean(X);
k=ones(1,sz(2));

for m=1:sz(2)
    x=X(:,m);
xi4=(1/n)*sum((x-xbar(m)).^4);
xi2=(1/n)*sum((x-xbar(m)).^2);
k(m)=(xi4/(xi2.^2));
end

if ~flag
    k=(((n-1)/((n-2)*(n-3)))*(((n+1)*k)-(3*(n-1))))+3;
end

end