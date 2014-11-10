function flag = isvec( v )
%
% flag = isvec( v )
%
% returns 1 if min(size(v))==1 & length(v)>1
%   and   0 otherwise
%
% by G.F. Margrave, October 1993
%
	if( min(size(v))==1 & length(v) > 1 )
		flag = 1;
	else
		flag = 0;
	end
