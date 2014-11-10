function X = jddiff(X,nd)
%DIFF	Differentiate or difference.
%	DIFF may be called numeric arguments.
%
%	For a numeric vector argument, DIFF computes differences.
%	DIFF(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)].
%	DIFF(X), for a matrix X, is the matrix of column differences,
%	   [X(2:n,:) - X(1:n-1,:)].
%	DIFF(X,nd) is the index seperation for the difference function.
%	For example if nd=5, DIFF(X) is [X(6)-X(1)  X(7)-X(2) ... X(n)-X(n-5)].
[m,n] = size(X);
if m == 1
   X = X(nd+1:n) - X(1:n-nd);
else
   X = X(nd+1:m,:) - X(1:m-nd,:);
end
