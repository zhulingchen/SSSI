% MAINFWIFREQCPMLFOR2DAW simulates the full waveform inversion (FWI) with
% 2-d acoustic wave in frequency domain based on the CPML absorbing
% boundary condition
%
% The FWI in frequency domain is used to solve the following problem: Given
% a smooth but obscure velocity model as a starting point and then minimize
% the least-squares misfit function defined by the differences at the
% receiver positions between the recorded seismic data and the modeled
% seismic data for each source-receiver pair of the seismic survey.
%
%
% System background
% ====================================================================================================
%
% The true velocity model m is approached by iteratively running the FWI on
% the currently estimated velocity model m'
%
% PDE for true field u_obs:                    (m(x)(d^2/dt^2) - Laplacian)u_obs(x, t; xs) = f(x, t; xs)
% PDE for calculated (modeled) field u_cal:    (m'(x)(d^2/dt^2) - Laplacian)u_cal(x, t; xs) = f(x, t; xs)
% u_obs(x, t; xs) is the observed (recorded) data at position x caused by
% the shot at position xs
% u_cal(x, t; xs) is the calculated (modeled) data at position x caused by
% the shot at position xs
%
% In frequency domain
% PDE for true field u_obs:                    (- m(x)w^2 - Laplacian)U_obs(x, jw; xs) = F(x, jw; xs)
% PDE for calculated (modeled) field u_cal:    (- m'(x)w^2 - Laplacian)U_cal(x, jw; xs) = F(x, jw; xs)
% whose solution is
% U_cal(x, jw; xs) = G'(x, jw; xs) * F(x, jw; xs)
% or
% U_cal(x, jw; xs) = A(jw; xs) \ (-F(x, jw; xs))
% where G'(y, jw; x) is the Green's function of m' from source xs to
% receiver x and A(jw; xs) is a discretization matrix (Helmholtz operator)
% in frequency domain that maps 2D stencils into its columns
%
%
% The cost function is:
% J = 1/2 * \sum_w \sum_xs \sum_xr |U_cal(xr, jw; xs) - U_obs(xr, jw; xs)|^2
%
% ====================================================================================================
%
%
% Purpose
% ====================================================================================================
%
% To find an estimate m' of m such that J is minimized
%
% The minimum of the misfit function J is sought in the vicinity of the
% starting model m_0', m_{k+1}' = m_k' + a_k * g_k where g_k is the
% gradient of J for m_k'
%
% ====================================================================================================
%
%
% Method
% ====================================================================================================
%
% J is minimized using quasi-Newton method
%
% ====================================================================================================
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;