close all;
clear;
clc;

% interface for testMatrixMultiply_mpiAgent
A = randn(256, 512);
B = randn(512, 256);
C = A * B;
C_mpi = testMatrixMultiply_mpiAgent(A, B);
[m, n] = size(C_mpi);
fprintf('m = %d, n = %d\n', m, n);
delta = C - C_mpi;
max(abs(delta(:)))

quit;