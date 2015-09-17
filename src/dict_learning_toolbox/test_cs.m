close all;
clear;
clc;

N = 128;
M = 32;

% Compressive Sampling
theta = 2*pi * rand(N, 1);
fftMat = 1/sqrt(N) * fft(eye(N, N));
csMat = diag(exp(1j*theta));
randSign = binornd(1, 0.5, [N, N]);
randSign(randSign == 0) = -1;
csMat = randSign .* csMat;
randSrcGrid = randperm(N);
randSrcGrid = sort(randSrcGrid(1:M));
csMat = csMat(randSrcGrid, :);

% test
tmpFreq = zeros(N, 1);
tmpFreq(20:20:100) = 1;
ifftMat = fftMat';
tmpTime = sqrt(N) * ifft(tmpFreq, N);

% BP (time domain)
tmpFreq_rec_bpdn = spg_bp(csMat * ifftMat, csMat * tmpTime);
figure; plot(real(tmpFreq_rec_bpdn));
% LASSO (time domain)
tmpFreq_rec_lasso = spg_lasso(csMat * ifftMat, csMat * tmpTime, norm(tmpFreq, 1));
figure; plot(real(tmpFreq_rec_lasso));

% BP (freq domain)
csMat = 1/sqrt(2) * (randn(M, N) + 1j * randn(M, N));
tmpFreq_rec_bpdn = spg_bp(csMat, csMat * tmpFreq);
figure; plot(real(tmpFreq_rec_bpdn));
% LASSO (freq domain)
tmpFreq_rec_lasso = spg_lasso(csMat, csMat * tmpFreq, norm(tmpFreq, 1));
figure; plot(real(tmpFreq_rec_lasso));

% PQN
% Set up Objective Function
csMat = randn(M, N);
idctMat = idct(eye(N, N));
tmpTime = idct(tmpFreq, N);
funObj = @(w)SquaredError(w, csMat * idctMat, csMat * tmpTime);
% Set up L1-Ball Projection
tau = norm(tmpFreq, 1);
funProj = @(w)sign(w).*projectRandom2C(abs(w), tau);
options.optTol = 1e-10;

fprintf('\nComputing optimal Lasso parameters...\n');
tmpFreq_rec_pqn = zeros(N, 1);
[tmpFreq_rec_pqn, r] = minConF_PQN_new(funObj, tmpFreq_rec_pqn, funProj, options);
figure; plot(real(tmpFreq_rec_pqn));