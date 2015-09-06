close all;
clear;
clc;

N = 100;
M = 20;
% Compressive Sampling
theta = 2*pi * rand(N, 1);
fftMat = 1/sqrt(N) * fft(eye(N, N));
csMat = fftMat' * diag(exp(1j*theta)) * fftMat;
randSign = binornd(1, 0.5, [N, N]);
randSign(randSign == 0) = -1;
csMat = randSign .* csMat;
randSrcGrid = randperm(N);
randSrcGrid = sort(randSrcGrid(1:M));
csMat = csMat(randSrcGrid, :);

% test
tmpFreq = zeros(N, 1);
tmpFreq(20:20:80) = 1;
idctMat = idct(eye(N, N));
tmpTime = idct(tmpFreq, N);

% BP (time domain)
tmpFreq_rec_bpdn = spg_bp(csMat * idctMat, csMat * tmpTime);
figure; plot(real(tmpFreq_rec_bpdn));
% LASSO (time domain)
tmpFreq_rec_lasso = spg_lasso(csMat * idctMat, csMat * tmpTime, norm(tmpFreq, 1));
figure; plot(real(tmpFreq_rec_lasso));

% BP (freq domain)
tmpFreq_rec_bpdn = spg_bp(csMat, csMat * tmpFreq);
figure; plot(real(tmpFreq_rec_bpdn));
% LASSO (freq domain)
tmpFreq_rec_lasso = spg_lasso(csMat, csMat * tmpFreq, norm(tmpFreq, 1));
figure; plot(real(tmpFreq_rec_lasso));

% PQN
% Set up Objective Function
csMat = randn(20, 100);
funObj = @(w)SquaredError(w,csMat * idctMat, csMat * tmpTime);
% Set up L1-Ball Projection
tau = 4;
funProj = @(w)sign(w).*projectRandom2C(abs(w),tau);

fprintf('\nComputing optimal Lasso parameters...\n');
wL1 = zeros(N, 1);
[wL1, r] = minConF_PQN_new(funObj,wL1,funProj);
figure; plot(real(wL1));