close all;
clear;
clc;

MASTER = 0;

% interface for testMatrixMultiply_mpiAgent
A = randn(256, 512);
B = randn(512, 256);
C = A * B;
[C_mpi, taskId] = testMatrixMultiply_mpiAgent(A, B);
if (taskId == MASTER)
    [m, n] = size(C_mpi);
    fprintf('(In Matlab interface) Master: m = %d, n = %d, taskId = %d\n', m, n, taskId);
    delta = C - C_mpi;
    max(abs(delta(:)))
end

quit;