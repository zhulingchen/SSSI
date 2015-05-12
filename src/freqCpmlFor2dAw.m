function [A, snapshot] = freqCpmlFor2dAw(model, source, w, nDiffOrder, nBoundary, dz, dx)
%
% FREQCPMLFOR2DAW solves the following equation in frequency domain
%
% -(w^2)/(v^2)*U(z, x, jw) - (d^2)U(z, x, jw)/dz^2 - (d^2)U(z, x, jw)/dx^2 = S(z, x, jw)
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = -S(z, x, jw)
% m*(w^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = -S(z, x, jw)
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)*F(jw)
%                                           |
%                                           V
%                                     A * U = -S
% for U(z, x, jw) with with Nonsplit Convolutional-PML (CPML) Absorbing Boundary Conditions
%
% input arguments
% model             velocity model (squared slowness)
% source            source vector (e.g., shots) in frequency domain
% w                 analog angular frequency \omega = [-pi, pi)/dt
% nDiffOrder        number of approximation order for differentiator operator
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
%
% output arguments
% A                 discretization matrix in frequency domain
% snapshot          pressure field u(z, x, jw) in frequency domain
%
% Reference:
% D. Komatitsch and R. Martin, An unsplit convolutional perfectly matched
% layer improved at grazing incidence for the seismic wave equation,
% Geophysics, Vol. 72 No. 5, pp. SM155-SM167, 2007
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

[nz, nx] = size(model);
velocity = sqrt(1./model);
nLength = nz * nx;
s = reshape(source, nLength, []);
coeff = dCoef(nDiffOrder, 's');
k = 2 * nDiffOrder - 1;


%% Absorbing boundary condition (ABC): Nonsplit Convolutional-PML (CPML)
ixb = 1:nBoundary;          % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft = dampPml(repmat(fliplr(ixb) * dx, nz, 1), velocity(:, ixb), nBoundary * dx);
xDampRight = dampPml(repmat(ixb * dx, nz, 1), velocity(:, ixb2), nBoundary * dx);
xDamp = [xDampLeft, zeros(nz, nx-2*nBoundary), xDampRight];

zDampDown = dampPml(repmat(ixb.' * dz, 1, nx), velocity(izb2, :), nBoundary * dz);
zDamp = [zeros(nz-nBoundary, nx); zDampDown];


%% fill the elements into A
% with boundary padding (fastest and most readable)
modelExt = padarray(model, [k, k], 'replicate');
xDampExt = padarray(xDamp, [k, k], 'replicate');
zDampExt = padarray(zDamp, [k, k], 'replicate');
% interior domain
[ix, iz] = meshgrid((1+k):(nx+k), (1+k):(nz+k));
idxRowInternal = (ix-1)*(nz+2*k) + iz;
idxColInternal = zeros(nLength * (4*k+1), 1);
valInternal = zeros(nLength * (4*k+1), 1);
% summed coefficients for each offset
c = zeros(1, 2*k+1);
for ii = 1:nDiffOrder
    for jj = 1:nDiffOrder
        iOffset1 = ii + jj - 1;
        iOffset2 = ii - jj;
        c(iOffset1+(k+1)) = c(iOffset1+(k+1)) + coeff(ii) * coeff(jj);
        c(iOffset2+(k+1)) = c(iOffset2+(k+1)) - coeff(ii) * coeff(jj);
        c(-iOffset1+(k+1)) = c(-iOffset1+(k+1)) + coeff(ii) * coeff(jj);
        c(-iOffset2+(k+1)) = c(-iOffset2+(k+1)) - coeff(ii) * coeff(jj);
    end
end
% serialized index of left, right, up, down neighbor elements and the
% center element itself
ii = 1;
for iOffset = [-k:-1, 1:k]
    % index of left, right neighbor elements
    idxColInternal((ii - 1) * nLength + (1:nLength)) = (ix-1+iOffset)*(nz+2*k) + iz;
    % values of left, right neighbor elements
    valInternal((ii - 1) * nLength + (1:nLength)) = c(iOffset+(k+1)) * (1j*w./(dx*(1j*w+xDampExt(idxRowInternal)))).^2;
    % index of up, down neighbor elements
    idxColInternal((2 * k + ii) * nLength + (1:nLength)) = (ix-1)*(nz+2*k) + (iz+iOffset);
    % values of up, down neighbor elements
    valInternal((2 * k + ii) * nLength + (1:nLength)) = c(iOffset+(k+1)) * (1j*w./(dz*(1j*w+zDampExt(idxRowInternal)))).^2;
    ii = ii + 1;
end
% index of the center element itself
idxColInternal((ii - 1) * nLength + (1:nLength)) = (ix-1)*(nz+2*k) + iz;
% value of the center element itself
valInternal((ii - 1) * nLength + (1:nLength)) = (modelExt(idxRowInternal) .* w^2) ...
    + c(k+1) * (1j*w./(dz*(1j*w+zDampExt(idxRowInternal)))).^2 ...
    + c(k+1) * (1j*w./(dx*(1j*w+xDampExt(idxRowInternal)))).^2;
% create sparse matrix
A = sparse(repmat(idxRowInternal(:), 4*k+1, 1), idxColInternal, valInternal);
% remove padding
A = A(idxRowInternal, idxRowInternal);


% % using double-nested for loops (much slower, cannot be used for higher-order approximation of staggered-grid finite difference)
% A2 = sparse(nLength, nLength);
% for iz = 1:nz
%     for ix = 1:nx
%         M = zeros(nz, nx);
%         if (iz-1 >= 1)
%             M(iz-1, ix) = (1j*w/(dz*(1j*w+zDamp(iz, ix))))^2; % up
%         end
%         if (ix-1 >= 1)
%             M(iz, ix-1) = (1j*w/(dx*(1j*w+xDamp(iz, ix))))^2; % left
%         end
%         M(iz, ix) = (model(iz, ix) * w^2) - 2*(1j*w/(dz*(1j*w+zDamp(iz, ix))))^2 - 2*(1j*w/(dx*(1j*w+xDamp(iz, ix))))^2; % center
%         if (ix+1 <= nx)
%             M(iz, ix+1) = (1j*w/(dx*(1j*w+xDamp(iz, ix))))^2; % right
%         end
%         if (iz+1 <= nz)
%             M(iz+1, ix) = (1j*w/(dz*(1j*w+zDamp(iz, ix))))^2; % down
%         end
%         row = reshape(M, 1, nLength);
%         A2((ix-1)*nz + iz, :) = row;
%     end
% end


%% A * U = -S, solve U(z, x, jw)

% tic;
% snapshot = umfpack2 (A, '\', (-s));
% toc;

% tic;
snapshot = A \ (-s);
% toc;

snapshot = reshape(snapshot, size(source));
