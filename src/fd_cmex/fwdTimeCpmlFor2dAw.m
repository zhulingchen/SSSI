function [data, snapshot] = fwdTimeCpmlFor2dAw(v, source, nDiffOrder, nBoundary, dz, dx, dt)
%
% FWDTIMECPMLFOR2DAW Simulate 2-d acoustic wave forward propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = zP + xP
% Update xPhi, zPhi, xA, zA, xPsi, zPsi, xP, zP and solve u(z, x, t) with Nonsplit Convolutional-PML (CPML)
%
% input arguments
% v(nz,nx)          velocity model
% source(nz,nx,nt)  source vector (e.g., shots)
% nDiffOrder        number of approximation order for differentiator operator
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
%
% output arguments
% data              received 2-d x - time signal on the surface
% snapshot          pressure field u(z, x, t) in time domain
%
% Reference: 
% D. Komatitsch and R. Martin, An unsplit convolutional perfectly matched
% layer improved at grazing incidence for the seismic wave equation,
% Geophysics, Vol. 72 No. 5, pp. SM155-SM167, 2007
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% It is a Matlab interface of its background C/Mex function to replace its
% Matlab function.
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

[data, snapshot] = fwdTimeCpmlFor2dAw_mex(v, source, nDiffOrder, nBoundary, dz, dx, dt);