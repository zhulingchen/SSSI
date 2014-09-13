function d = dampPml(u, v, L)
% DAMPPML Generate the model for damping parameter
% u = x or z, representing the distance between current position (in PML)
% and PML inner boundary
% 
% v, acoustic wave velocity
% L, the PML thickness
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

R = 1e-6; % a constant chosen as 1e-6 to 1e-3
d0 = -(3 * v)/(2 * L) * log(R);
d = d0 .* (u / L).^2;



