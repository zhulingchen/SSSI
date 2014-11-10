function polarization = getpol(x)
% This function will get the degree of polarization from a 3C time sequence
% by p = ((l1-l2)^2+(l1-l3)^2+(l2-l3)^2)/(2*(l1+l2+l3))
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lijun Zhu (gatechzhu@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


% check dimension of x
x = squeeze(x);
if size(x,1) ~= 3
    x = x';
    if size(x,1)~=3
        error('Wrong dimension');
    end
end
    
eigenvalue = eig(x*x');
l1 = eigenvalue(1);
l2 = eigenvalue(2);
l3 = eigenvalue(3);
polarization = ((l1-l2)^2+(l1-l3)^2+(l2-l3)^2)/(2*(l1+l2+l3)^2);
end
