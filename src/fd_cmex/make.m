close all;
clear;
clc;

fprintf('Compiling...\n');
if (isunix) % Linux / MacOS
    mex COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3" diffOperator_mex.c
    mex COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3" fwdTimeCpmlFor2dAw_mex.c finiteDifference.c
    mex COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3" rvsTimeCpmlFor2dAw_mex.c finiteDifference.c
else        % Windows
    mex diffOperator_mex.c
    mex fwdTimeCpmlFor2dAw_mex.c finiteDifference.c
    mex rvsTimeCpmlFor2dAw_mex.c finiteDifference.c
end

fprintf('Compiling with OpenMP...\n');
if (isunix) % Linux / MacOS
    mex COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" fwdTimeCpmlFor2dAw_openmp_mex.c finiteDifference.c
    mex COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" rvsTimeCpmlFor2dAw_openmp_mex.c finiteDifference.c
else        % Windows
    mex CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" fwdTimeCpmlFor2dAw_openmp_mex.c finiteDifference.c
    mex CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" rvsTimeCpmlFor2dAw_openmp_mex.c finiteDifference.c
end

if (isunix) % Linux / MacOS
    fprintf('Compiling with OpenMPI...\n');
    mex CC='/usr/local/bin/mpicc' LD='/usr/local/bin/mpicc' CFLAGS='\$CFLAGS -std=c99' fwdTimeCpmlFor2dAw_openmpi_mex.c finiteDifference.c
end
fprintf('Compiling complete!\n');
