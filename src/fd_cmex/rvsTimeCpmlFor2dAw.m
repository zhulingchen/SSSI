function [model, snapshot] = rvsTimeCpmlFor2dAw(v, data, nDiffOrder, nBoundary, dz, dx, dt)
%
% RVSTIMECPMLFOR2DAW Simulate 2-d acoustic wave reverse propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = zP + xP
% Update xPhi, zPhi, xA, zA, xPsi, zPsi, xP, zP and solve u(z, x, t) with Nonsplit Convolutional-PML (CPML)
%
% input arguments
% v(nz,nx)          velocity model
% data(nx,nt)       received data on the surface
% nDiffOrder        number of approximation order for differentiator operator
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
%
% output arguments
% model             final pressure wavefield u(z, x, 1) in reverse time domain
% snapshot          pressure field u(z, x, t) in time domain
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% It is a Matlab interface of its background C/Mex function to replace its
% Matlab function.
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

[model, snapshot] = rvsTimeCpmlFor2dAw_mex(v, data, nDiffOrder, nBoundary, dz, dx, dt);