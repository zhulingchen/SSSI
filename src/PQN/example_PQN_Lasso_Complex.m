function example_PQN_Lasso_Complex
% solve min_x ||R*fft(x)-b||^2 s.t. ||x||_1 <= tau

close all;
clear;
clc;

addpath(genpath('./'));

m = 128;
n = 512;
R = randn(m, n) + 1j * randn(m, n);
R = R / sqrt(m);
x = 10 * randn(n,1).*(rand(n,1) > 0.9);

X = 1/sqrt(n) * fft(x, n);

F = 1/sqrt(n) * fft(eye(n, n));

b = R * X;
tau = norm(x, 1);

%% SPGL1
xL1_spgLasso = spg_lasso(R * F, b, tau);

%% minConF
options.verbose = 2;
options.optTol = 1e-8;
options.SPGoptTol = 1e-25;
options.adjustStep = 1;
options.bbInit = 1;
options.maxIter = 2000;

func = @(zRealImag) misfitFunc([real(R) -imag(R); imag(R) real(R)], zRealImag, [real(b); imag(b)]);

% L2-norm minimization with L1-norm regularization, i.e., LASSO
tau = norm(x, 1);
funProj = @(zRealImag) complexProject(zRealImag, tau);

% % L2-norm minimization without L1-norm regularization
% funProj = @(x) boundProject(x, [-inf(2*n, 1)], [inf(2*n, 1)]);

xL1_pqnRealImag = minConF_PQN_new(func, [zeros(n, 1); zeros(n, 1)], funProj, options);
xL1_pqn = xL1_pqnRealImag(1:n) + 1j*xL1_pqnRealImag(n+1:2*n);

end

function [value, grad] = misfitFunc(R, x, b)

n = length(x);

Xreal = 1/sqrt(n/2) * ( real(fft(x(1:n/2), n/2)) - imag(fft(x(n/2+1:n), n/2)) );
Ximag = 1/sqrt(n/2) * ( imag(fft(x(1:n/2), n/2)) + real(fft(x(n/2+1:n), n/2)) );
X = [Xreal; Ximag];

value = 1/2 * sum((R*X-b).^2);

grad = R'*(R*X-b);

gradReal = sqrt(n/2) * ( real(ifft(grad(1:n/2), n/2)) - imag(ifft(grad(n/2+1:n), n/2)) );
gradImag = sqrt(n/2) * ( imag(ifft(grad(1:n/2), n/2)) + real(ifft(grad(n/2+1:n), n/2)) );

grad = [gradReal; gradImag];

end