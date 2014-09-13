function [A, snapshot] = freqSolver(v, source, w, nBoundary, dz, dx)
%
% FREQSOLVER solves the following equation in frequency domain
% 
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = S(z, x, jw)
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)*F(jw)
%                                           |
%                                           V
%                                     A * U = S
% for U(z, x, jw) with Enguist-Majda Absorbing Boundary Conditions
%
% input arguments
% v(nz,nx)          velocity model
% source(nz,nx)     source vector (e.g., shots)
% w                 analog angular frequency \omega = [-pi, pi)/dt
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% 
% output arguments
% A                 discretization matrix in frequency domain
% snapshot          pressure field u(z, x, jw) in frequency domain
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
velocityModel = v(1:end-nBoundary, nBoundary+1:end-nBoundary);
[nz, nx] = size(velocityModel);
nLength = nz * nx;
s = source(1:end-nBoundary, nBoundary+1:end-nBoundary);
s = reshape(s, nLength, 1);

A = sparse(nLength, nLength);


%% fill the elements into A
% using the method introduced in the paper "Frequency-Domain Modeling
% Techniques for the Scalar Wave Equation : An Introduction"
% Author: Ajo-Franklin, Jonathan B.

% domain interior
for iz = 2:nz-1
    for ix = 2:nx-1
        M = zeros(nz, nx);
        M(iz-1, ix) = 1/dz^2;
        M(iz, ix-1) = 1/dx^2;
        M(iz, ix) = (w^2/velocityModel(iz, ix)^2) - 2*(1/dz^2 + 1/dx^2);
        M(iz, ix+1) = 1/dx^2;
        M(iz+1, ix) = 1/dz^2;
        row = reshape(M, 1, nLength);
        A((ix-1)*nz + iz, :) = row;
    end
end

% Second-order Enguist-Majda Absorbing Boundary Conditions
% see article "absorbing boundary conditions - University of Kentucky"
% http://www.engr.uky.edu/~gedney/courses/ee624/notes/EE624_Notes6.pdf

% top boundary (iz = 1)
for ix = 2:nx-1
    M = zeros(nz, nx);
    M(1, ix-1) = velocityModel(1, ix) / (2 * dx^2);
    M(1, ix) = -1j*w/dz + w^2/velocityModel(1, ix) - velocityModel(1, ix)/dx^2;
    M(1, ix+1) = velocityModel(1, ix) / (2 * dx^2);
    M(2, ix) = 1j*w/dz;
    row = reshape(M, 1, nLength);
    A((ix-1)*nz + 1, :) = row;
end

% bottom boundary (iz = nz)
for ix = 2:nx-1
    M = zeros(nz, nx);
    M(nz-1, ix) = 0;
    M(nz, ix-1) = -velocityModel(nz, ix) / (2 * dx^2);
    M(nz, ix) = -1j*w/dz - w^2/velocityModel(nz, ix) + velocityModel(nz, ix)/dx^2;
    M(nz, ix+1) = -velocityModel(nz, ix) / (2 * dx^2);
    row = reshape(M, 1, nLength);
    A((ix-1)*nz + nz, :) = row;
end

% left boundary (ix = 1)
for iz = 2:nz-1
    M = zeros(nz, nx);
    M(iz-1, 1) = velocityModel(iz, 1) / (2 * dz^2);
    M(iz, 1) = -1j*w/dx + w^2/velocityModel(iz, 1) - velocityModel(iz, 1)/dz^2;
    M(iz, 2) = 1j*w/dx;
    M(iz+1, 1) = velocityModel(iz, 1) / (2 * dz^2);
    row = reshape(M, 1, nLength);
    A(iz, :) = row;
end

% right boundary (ix = nx)
for iz = 2:nz-1
    M = zeros(nz, nx);
    M(iz-1, nx) = -velocityModel(iz, nx) / (2 * dz^2);
    M(iz, nx-1) = 0;
    M(iz, nx) = -1j*w/dx - w^2/velocityModel(iz, nx) + velocityModel(iz, nx)/dz^2;
    M(iz+1, nx) = -velocityModel(iz, nx) / (2 * dz^2);
    row = reshape(M, 1, nLength);
    A((nx-1)*nz + iz, :) = row;
end

% top left corner (iz = 1, ix = 1)
M = zeros(nz, nx);
M(1, 1) = -1j*w/dz - 1j*w/dx + 2*w^2/velocityModel(1, 1) - velocityModel(1, 1)/dx^2 - velocityModel(1, 1)/dz^2;
M(1, 2) = 1j*w/dx + velocityModel(1, 1) / (2 * dx^2);
M(2, 1) = 1j*w/dz + velocityModel(1, 1) / (2 * dz^2);
row = reshape(M, 1, nLength);
A(1, :) = row;

% top right corner (iz = 1, ix = nx)
M = zeros(nz, nx);
M(1, nx-1) = velocityModel(1, nx) / (2 * dx^2);
M(1, nx) = -1j*w/dz - 1j*w/dx - velocityModel(1, nx)/dx^2 + velocityModel(1, nx)/dz^2;
M(2, nx) = 1j*w/dz - velocityModel(1, nx) / (2 * dz^2);
row = reshape(M, 1, nLength);
A((nx-1)*nz + 1, :) = row;

% bottom left corner (iz = nz, ix = 1)
M = zeros(nz, nx);
M(nz-1, 1) = velocityModel(nz, 1) / (2 * dz^2);
M(nz, 1) = -1j*w/dz - 1j*w/dx + velocityModel(nz, 1)/dx^2 - velocityModel(nz, 1)/dz^2;
M(nz, 2) = 1j*w/dx - velocityModel(nz, 1) / (2 * dx^2);
row = reshape(M, 1, nLength);
A(nz, :) = row;

% bottom right corner (iz = nz, ix = nx)
M = zeros(nz, nx);
M(nz-1, nx) = -velocityModel(nz, nx) / (2 * dz^2);
M(nz, nx-1) = -velocityModel(nz, nx) / (2 * dx^2);
M(nz, nx) = -1j*w/dz - 1j*w/dx - 2*w^2/velocityModel(nz, nx) + velocityModel(nz, nx)/dx^2 + velocityModel(nz, nx)/dz^2;
row = reshape(M, 1, nLength);
A((nx-1)*nz + nz, :) = row;

% % First-order Enguist-Majda Absorbing Boundary Conditions
% 
% % top boundary (iz = 1)
% for ix = 2:nx-1
%     M = zeros(nz, nx);
%     M(1, ix-1) = 0;
%     M(1, ix) = -1/dz - 1j*w/velocityModel(1, ix);
%     M(1, ix+1) = 0;
%     M(2, ix) = 1/dz;
%     row = reshape(M, 1, nLength);
%     A((ix-1)*nz + 1, :) = row;
% end
% 
% % bottom boundary (iz = nz)
% for ix = 2:nx-1
%     M = zeros(nz, nx);
%     M(nz-1, ix) = 1/dz;
%     M(nz, ix-1) = 0;
%     M(nz, ix) = -1/dz - 1j*w/velocityModel(nz, ix);
%     M(nz, ix+1) = 0;
%     row = reshape(M, 1, nLength);
%     A((ix-1)*nz + nz, :) = row;
% end
% 
% % left boundary (ix = 1)
% for iz = 2:nz-1
%     M = zeros(nz, nx);
%     M(iz-1, 1) = 0;
%     M(iz, 1) = -1/dx - 1j*w/velocityModel(iz, 1);
%     M(iz, 2) = 1/dx;
%     M(iz+1, 1) = 0;
%     row = reshape(M, 1, nLength);
%     A(iz, :) = row;
% end
% 
% % right boundary (ix = nx)
% for iz = 2:nz-1
%     M = zeros(nz, nx);
%     M(iz-1, nx) = 0;
%     M(iz, nx-1) = 1/dx;
%     M(iz, nx) = -1/dx - 1j*w/velocityModel(iz, nx);
%     M(iz+1, nx) = 0;
%     row = reshape(M, 1, nLength);
%     A((nx-1)*nz + iz, :) = row;
% end
% 
% % dp/dn at corner locations are the sum of the two orthogonal boundary stencils.
% 
% % top left corner (iz = 1, ix = 1)
% M = zeros(nz, nx);
% M(1, 1) = -1/dz - 1/dx - 2j*w/velocityModel(1, 1);
% M(1, 2) = 1/dx;
% M(2, 1) = 1/dz;
% row = reshape(M, 1, nLength);
% A(1, :) = row;
% 
% % top right corner (iz = 1, ix = nx)
% M = zeros(nz, nx);
% M(1, nx-1) = 1/dx;
% M(1, nx) = -1/dz - 1/dx - 2j*w/velocityModel(1, nx);
% M(2, nx) = 1/dz;
% row = reshape(M, 1, nLength);
% A((nx-1)*nz + 1, :) = row;
% 
% % bottom left corner (iz = nz, ix = 1)
% M = zeros(nz, nx);
% M(nz-1, 1) = 1/dz;
% M(nz, 1) = -1/dz - 1/dx - 2j*w/velocityModel(nz, 1);
% M(nz, 2) = 1/dx;
% row = reshape(M, 1, nLength);
% A(nz, :) = row;
% 
% % bottom right corner (iz = nz, ix = nx)
% M = zeros(nz, nx);
% M(nz-1, nx) = 1/dz;
% M(nz, nx-1) = 1/dx;
% M(nz, nx) = -1/dz - 1/dx - 2j*w/velocityModel(nz, nx);
% row = reshape(M, 1, nLength);
% A((nx-1)*nz + nz, :) = row;
%
% % dp/dn at corner locations are at 45 degrees from either grid axis, in a rotated coordinate frame
%
% % top left corner (iz = 1, ix = 1)
% M = zeros(nz, nx);
% M(1, 1) = -1/sqrt(dz^2 + dx^2) - 1j*w/velocityModel(1, 1);
% M(2, 2) = 1/sqrt(dz^2 + dx^2);
% row = reshape(M, 1, nLength);
% A(1, :) = row;
% 
% % top right corner (iz = 1, ix = nx)
% M = zeros(nz, nx);
% M(1, nx) = -1/sqrt(dz^2 + dx^2) - 1j*w/velocityModel(1, nx);
% M(2, nx-1) = 1/sqrt(dz^2 + dx^2);
% row = reshape(M, 1, nLength);
% A((nx-1)*nz + 1, :) = row;
% 
% % bottom left corner (iz = nz, ix = 1)
% M = zeros(nz, nx);
% M(nz-1, 2) = 1/sqrt(dz^2 + dx^2);
% M(nz, 1) = -1/sqrt(dz^2 + dx^2) - 1j*w/velocityModel(nz, 1);
% row = reshape(M, 1, nLength);
% A(nz, :) = row;
% 
% % bottom right corner (iz = nz, ix = nx)
% M = zeros(nz, nx);
% M(nz-1, nx-1) = 1/sqrt(dz^2 + dx^2);
% M(nz, nx) = -1/sqrt(dz^2 + dx^2) - 1j*w/velocityModel(nz, nx);
% row = reshape(M, 1, nLength);
% A((nx-1)*nz + nz, :) = row;

%% A * U = S, solve U(z, x, jw)

% snapshot = umfpack2 (A, '\', s);

snapshot = A \ s;

% snapshot = reshape(snapshot, nz, nx);
