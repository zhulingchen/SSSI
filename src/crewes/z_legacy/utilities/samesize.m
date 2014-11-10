function flag = samesize(a,b)
%
% flag = samesize(a,b)
%
% returns 1 if a and b are exactly the same size and zero otherwise
%
	flag = prod(size(a)==size(b));
