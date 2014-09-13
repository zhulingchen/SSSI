function c = dCoef(o, s)
% DCOEF is used to calculate differential coefficient.
% o is the half precision order of difference.
% s stands for difference type. It have two values,'r' and 's'.
% 'r' indicates regular-grid differnece;'s' indicate staggered-grid difference
%
% Author: Yan Xinfei(fulcrum)
% Finishing Time: 2007/10/20
%
% Usage:
%       c = DCoef(o) ruturn 2*o order coefficient of regular-grid difference.
%       c = DCoef(o,'r') the same as above
%       c = DCoef(o,'s') ruturn 2*o order coefficient of stragged-grid difference.

c = zeros(o, 1);
A = zeros(o);
if nargin < 1
    error('Please enter input arguments!');
elseif nargin == 1
    b = [1/2; zeros(o-1, 1)];
    for i=1:o
        for j=1:o
            A(i,j)=j^(2*i-1);
        end
    end
elseif nargin == 2
    if s=='r'
        b = [1/2; zeros(o-1, 1)];
        for i=1:o
            for j=1:o
                A(i,j)=j^(2*i-1);
            end
        end
    elseif s=='s'
        b = [1; zeros(o-1, 1)];
        for i=1:o
            for j=1:o
                A(i,j)=(2*j-1)^(2*i-1);
            end
        end
    else
        error('Input argument s can only be ''r'' or ''s''!');
    end
else
    error('Too many input aruguments!');
end

c = A \ b;