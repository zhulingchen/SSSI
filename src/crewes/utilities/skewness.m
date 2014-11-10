function s=skewness(X,flag)
% s=skewness(X), s=skewness(X,flag)
% Skewness is a measure of the asymmetry of the probablility distribution
% of a real random variable.  It is used in proabablitiy and statistics.
%
% This function will return the skewness of set of data.  If the data is a
% 2-D matrix it will return a value for each column.
%
% if flag=0 a scaling is applied. flag=1; is the default.

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
s=ones(1,sz(2));

for m=1:sz(2)
    x=X(:,m);
xi3=(1/n)*sum((x-xbar(m)).^3);
xi2=sqrt((1/n)*sum((x-xbar(m)).^2));
s(m)=xi3/(xi2.^3);
end

if ~flag
    s=sqrt((n*(n-1))/(n-2))*s;
end

end
