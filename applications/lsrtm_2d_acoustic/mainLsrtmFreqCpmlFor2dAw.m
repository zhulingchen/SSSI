% MAINLSRTMFREQCPMLFOR2DAW simulates the least-squares reverse time
% migration (LSRTM) with 2-d acoustic wave in frequency domain based on the
% CPML absorbing boundary condition, using Gauss-Newton method where the
% Hessian matrix approximated by its diagonal terms
%
% The LSRTM in frequency domain is used to solve the following problem:
% Given a smooth but obscure velocity model and the received data on
% surface, the true but unknown velocity model is to be approximated by
% estimating the scatter field during iterations.
%
% System background
% ====================================================================================================
%
% m = m_0 + epsilon * delta_m
% True field u: (m(x)(d^2/dt^2) - Laplacian)u(x, t; xs) = -f(x, t; xs)
% Incident field u_0: (m_0(x)(d^2/dt^2) - Laplacian)u_0(x, t; xs) = -f(x, t; xs)
% u = u_0 + u_sc, u_sc is scattered field
%
% Therefore,
% (m_0(x)(d^2/dt^2) - Laplacian)u_sc(y, t; xs) = -epsilon * delta_m(x) * (d^2/dt^2)u(x, t; xs)
% so that
% u_sc = -epsilon * G(\delta_m * (d^2/dt^2)u)
% u = (I + epsilon*A)^(-1)u_0
%   = u_0 - epsilon * A * u_0 + epsilon^2 * A^2 * u_0 - ...
%   = u_0 + epsilon * u_1 + epsilon^2 * u_2 + ...
%
% Therefore, u_1 = -A * u_0, where
% A(f) = G(\delta_m*(d^2/dt^2)f)
% and u_1 satisfies
%
% (m_0(x)(d^2/dt^2) - Laplacian)u_1(y, t; xs) = -delta_m(x) * (d^2/dt^2)u_0(x, t; xs) in time domain
% (-m_0(x)w^2 - Laplacian)U_1(y, jw; xs) = delta_m(x) * w^2 * U_0(x, jw; xs) in frequency domain
% and the solution of U_1 can be written in a linear form as
% U_1 = L * \delta_m
%
% The cost function is:
% J = 1/2 * \sum\w \sum\xs \sum\xr |U_1(xs, xr, jw) - \delta_D(xs, xr, jw)|^2 + lambda * |\delta_m(x)|^2
%
% ====================================================================================================
%
% Purpose
% ====================================================================================================
%
% To find an optimized delta_m(x) such that J is minimized and update the
% velocity model m
%
% ====================================================================================================
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;


%% Full Wave Inversion Example
% in frequency domain

ALPHA = 0.75;
DELTA = 1e-4;
EPSILON = 1;
FREQTHRES = 2;
MAXITER = 20;


%% Set path
run([fileparts(pwd), '/setpath']);


%% Read in velocity model data
filenameVelocityModel = [model_data_path, '/velocityModel.mat'];
[pathVelocityModel, nameVelocityModel] = fileparts(filenameVelocityModel);
load(filenameVelocityModel); % velocityModel
[nz, nx] = size(velocityModel);

% smooth velocity model used average filter
filenameVelocityModelSmooth = [model_data_path, '/velocityModelSmooth.mat'];
load(filenameVelocityModelSmooth); % velocityModelSmooth

nBoundary = 20;
% % a more smooth velocity model for FWI
% VS = extBoundary(velocityModelSmooth, nBoundary, 2);
% VS = [repmat(VS(1, :), nBoundary, 1); VS];
% nAvgSize = [1, 1];
% hImageSmooth = fspecial('average', nAvgSize);
% VS = imfilter(VS, hImageSmooth);
% velocityModelSmooth = VS(nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary);

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;

% grids and positions of shot array
shotArrType = 'uniform';
idxShotArrLeft = 1;
idxShotArrRight = nx;
nShots = nx;
if (strcmpi(shotArrType, 'uniform'))
    xShotGrid = (idxShotArrLeft:ceil((idxShotArrRight - idxShotArrLeft + 1)/nShots):idxShotArrRight);
elseif (strcmpi(shotArrType, 'random'))
    xShotGrid = (idxShotArrLeft:idxShotArrRight);
    xShotGrid = sort(xShotGrid(randperm(idxShotArrRight - idxShotArrLeft + 1, nShots)));
else
    error('Shot array deployment type error!');
end
xShot = xShotGrid * dx;

shotWatchList = [1, ceil(nShots/2), nShots];

% grids and positions of receiver array
recArrType = 'uniform';
idxRecArrLeft = 1;
idxRecArrRight = nx;
nRecs = nx;
if (strcmpi(recArrType, 'uniform'))
    xRecGrid = (idxRecArrLeft:ceil((idxRecArrRight - idxRecArrLeft + 1)/nRecs):idxRecArrRight);
elseif (strcmpi(recArrType, 'random'))
    xRecGrid = (idxRecArrLeft:idxRecArrRight);
    xRecGrid = sort(xRecGrid(randperm(idxRecArrRight - idxRecArrLeft + 1, nRecs)));
else
    error('Receiver array deployment type error!');
end
xRec = xRecGrid * dx;

xShotAndRecGrid = union(xShotGrid, xRecGrid);
nShotsAndRecs = length(xShotAndRecGrid);

% plot the velocity model
figure(1);
imagesc(x, z, velocityModel);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
colormap(seismic);


%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences with a continuous source function
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = ALPHA * (dz/vmax/sqrt(2));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;
nfft = 2^(nextpow2(nt));
dw = 2*pi/nfft;
w = (-pi:dw:pi-dw)/dt; % analog angular frequency \omega = [-pi, pi)/dt

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 2);
VS = extBoundary(velocityModelSmooth, nBoundary, 2);

% dimension of frequency-domain solution
nLength = nz * nx;
nLengthWithBoundary = (nz + nBoundary) * (nx + 2*nBoundary);

% number of approximation order for differentiator operator
nDiffOrder = 1;

% Define analog frequency parameter for ricker wavelet
f = 20;
% f = w(550)/(2*pi);


%% Shot data recording at the surface
% generate shot signal
rw1dTime = zeros(1, nt);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + ricker(f(ifreq), nt, dt);
end
rw1dFreq = fftshift(fft(rw1dTime, nfft), 2);
% find active frequency set with FFT amplitude larger than the threshold
activeW = find(abs(rw1dFreq) > FREQTHRES);
activeW = activeW(activeW > nfft / 2 + 1); % choose f > 0Hz

dataTrueFreq = zeros(nRecs, nShots, nfft);
dataDeltaFreq = zeros(nRecs, nShots, nfft);
% receiver positions on extended velocity model
xr = xRecGrid + nBoundary;

for ixs = 1:nShots %21:nx+20 % shot loop

    curXsPos = xShotGrid(ixs) + nBoundary; % shot position on x

    % generate shot signal
    sourceTime = zeros([size(V), nt]);
    sourceTime(1, curXsPos, :) = reshape(rw1dTime, 1, 1, nt);

    % generate shot record
    tic;
    [dataTrue, ~] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    [dataSmooth, ~] = fwdTimeCpmlFor2dAw(VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    timeForward = toc;
    fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, elapsed time = %fs\n', curXsPos-nBoundary, x(curXsPos-nBoundary), timeForward);

    dataTrue = dataTrue(xr, :);
    dataSmooth = dataSmooth(xr, :);

    dataTrueFreq(:, ixs, :) = fftshift(fft(dataTrue, nfft, 2), 2);
    dataDeltaFreq(:, ixs, :) = fftshift(fft(dataTrue - dataSmooth, nfft, 2), 2);

end % end shot loop

% save received surface data
filenameDataTrueFreq = [pathVelocityModel, '/dataTrueFreq.mat'];
save(filenameDataTrueFreq, 'dataTrueFreq', '-v7.3');

filenameDataDeltaFreq = [pathVelocityModel, '/dataDeltaFreq0.mat'];
save(filenameDataDeltaFreq, 'dataDeltaFreq', '-v7.3');

% clear variables and functions from memory
clear('dataTrueFreq');
clear('dataDeltaFreq');
clear('sourceTime');

%% Full wave inversion (FWI)
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + s(z, x, t) = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)F(jw)
%                                           |
%                                           V
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = S(z, x, jw)
%
% Green's function is the impulse response of the wave equation.

% generate impulse shot signal
impTime = [1, zeros(1, nt-1)];
impFreq = fftshift(fft(impTime, nfft), 2);

modelOld = zeros(nz, nx);
modelNew = 1./VS(1:end-nBoundary, nBoundary+1:end-nBoundary).^2;

% shot positions on extended velocity model
xs = xShotGrid + nBoundary;


%% Start a pool of Matlab workers
numCores = feature('numcores');
if isempty(gcp('nocreate')) % checking to see if my pool is already open
    myPool = parpool(numCores);
end

iter = 1;
while(norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro') > DELTA && iter <= MAXITER)
    
    modelOld = modelNew;
    vmOld = 1./sqrt(modelOld);
    vmOld = extBoundary(vmOld, nBoundary, 2);
    load(filenameDataDeltaFreq);
    
    % an approximated diagonal Hessian matrix
    hessianDiag = zeros(nLength, 1);
    % hessianMat = zeros(nLength, nLength);
    % migrated image
    mig = zeros(nLength, 1);
    
    figure(1);
    imagesc(x, z, vmOld(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Previous Velocity Model');
    colormap(seismic);
    
    % update the velocity model with least-squares
    parfor idx_w = 1:length(activeW)
        
        iw = activeW(idx_w);
        
        fprintf('Processing f(%d) = %fHz ... ', iw, w(iw)/(2*pi));
        tic;
        [A, ~] = freqCpmlFor2dAw(vmOld, zeros(size(V)), w(iw), nDiffOrder, nBoundary, dz, dx);
        
        % Green's function for every shot
        sourceFreq = zeros(nLengthWithBoundary, nShots);
        sourceFreq((xs-1)*(nz+nBoundary)+1, :) = impFreq(iw) * eye(nShots, nShots);
        greenFreqForShot = A \ sourceFreq;
        % remove external boundaries
        greenFreqForShot = reshape(greenFreqForShot, nz + nBoundary, nx + 2*nBoundary, nShots);
        greenFreqForShot = greenFreqForShot(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
        greenFreqForShot = reshape(greenFreqForShot, nLength, nShots);
        
        % Green's function for every receiver
        sourceFreq = zeros(nLengthWithBoundary, nRecs);
        sourceFreq((xr-1)*(nz+nBoundary)+1, :) = impFreq(iw) * eye(nRecs, nRecs);
        greenFreqForRec = A \ sourceFreq;
        % remove external boundaries
        greenFreqForRec = reshape(greenFreqForRec, nz + nBoundary, nx + 2*nBoundary, nRecs);
        greenFreqForRec = greenFreqForRec(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
        greenFreqForRec = reshape(greenFreqForRec, nLength, nRecs);
        
        % calculate the pseudo-Hessian matrix, which is the diagonal elements of Hessian matrix
        hessianDiag = hessianDiag + w(iw)^4 * abs(rw1dFreq(iw))^2 ...
            * sum(abs(greenFreqForShot).^2 .* (abs(greenFreqForRec).^2 * ones(nRecs, nShots)), 2);
        
        % % the true Hessian matrix
        % hessianMat = hessianMat + w(iw)^4 * abs(rw1dFreq(iw))^2 ...
        %     * (greenFreqForShot * greenFreqForShot') .* (greenFreqForRec * greenFreqForRec');
        
        % for ixs = 1:nShots
        %     for ixr = 1:nRecs
        %         hessianDiag2 = hessianDiag2 + real(w(iw)^4 * abs(rw1dFreq(iw))^2 ...
        %             * abs(greenFreqForShot(:, ixs)).^2 .* abs(greenFreqForRec(:, ixr)).^2);
        %     end
        % end
        
        % calculate the migrated image using dataDeltaFreq
        mig = mig + w(iw)^2 * (-rw1dFreq(iw)) * sum(greenFreqForShot .* (greenFreqForRec * conj(dataDeltaFreq(:, :, iw))), 2);
        
        % for ixs = 1:nShots
        %     for ixr = 1:nRecs
        %         grad2 = grad2 + real(w(iw)^2 * (-rw1dFreq(iw)) * conj(dataDeltaFreq(ixr, ixs, iw)) ...
        %             * greenFreqForShot(:, ixs) .* greenFreqForRec(:, ixr));
        %     end
        % end
        
        timePerFreq = toc;
        fprintf('elapsed time = %fs\n', timePerFreq);
        
    end
    
    % save the pseudo-Hessian matrix
    filenameHessianDiag = [pathVelocityModel, sprintf('/hessianDiag%d.mat', iter)];
    save(filenameHessianDiag, 'hessianDiag', '-v7.3');
    
    % save the migrated image
    filenameMig = [pathVelocityModel, sprintf('/mig%d.mat', iter)];
    save(filenameMig, 'mig', '-v7.3');
    
    lambda = 5 * max(hessianDiag);
    
    % updated model
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + EPSILON * (real(mig) ./ real(hessianDiag + lambda * ones(nLength, 1)));
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    figure(2);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model');
    colormap(seismic);
    % save current updated velocity model
    filenameVmNew = [pathVelocityModel, sprintf('/vmNew%d.mat', iter)];
    save(filenameVmNew, 'vmNew', 'modelNew', '-v7.3');
    
    
    %% debug begin
    %     % this code fragment implemented the original matrix L and bigL, which
    %     % is the stack of L's with respect to all xr and xs's. The code and its
    %     % results are proved to be the same with the above code but runs much
    %     % slower than that due to lots of large matrix multiplications. You can
    %     % regard this code fragment as the direct implementation of the
    %     % gradient and Hessian matrix of the cost function. We leave them here
    %     % for reference.
    %
    %     load(filenameDataDeltaFreq);
    %     L = zeros(1, nLength);
    %     bigL = zeros(nShots * nRecs, nLength);
    %     hessianTrue = zeros(nLength, nLength);
    %     hessianTrue2 = zeros(nLength, nLength);
    %     mig = zeros(nLength, 1);
    %     mig2 = zeros(nLength, 1);
    %     u1 = zeros(nShots, nRecs);
    %     u1_bak = zeros(nShots, nRecs);
    %     gradVerify = zeros(nLength, 1);
    %     gradVerify_bak = zeros(nLength, 1);
    %     for iw = activeW
    %
    %         if (iw == nfft / 2 + 1)
    %             % skip f = 0Hz
    %             continue;
    %         end
    %
    %         fprintf('Processing f(%d) = %fHz ... ', iw, w(iw)/(2*pi));
    %         tic;
    %         [A, ~] = freqCpmlFor2dAw(vmOld, zeros(size(V)), w(iw), nDiffOrder, nBoundary, dz, dx);
    %
    %         % Green's function for every shot
    %         sourceFreq = zeros(nLengthWithBoundary, nShots);
    %         sourceFreq((xs-1)*(nz+nBoundary)+1, :) = impFreq(iw) * eye(nShots, nShots);
    %         greenFreqForShot = A \ sourceFreq;
    %         % remove external boundaries
    %         greenFreqForShot = reshape(greenFreqForShot, nz + nBoundary, nx + 2*nBoundary, nShots);
    %         greenFreqForShot = greenFreqForShot(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
    %         greenFreqForShot = reshape(greenFreqForShot, nLength, nShots);
    %
    %         % Green's function for every receiver
    %         sourceFreq = zeros(nLengthWithBoundary, nRecs);
    %         sourceFreq((xr-1)*(nz+nBoundary)+1, :) = impFreq(iw) * eye(nRecs, nRecs);
    %         greenFreqForRec = A \ sourceFreq;
    %         % remove external boundaries
    %         greenFreqForRec = reshape(greenFreqForRec, nz + nBoundary, nx + 2*nBoundary, nRecs);
    %         greenFreqForRec = greenFreqForRec(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
    %         greenFreqForRec = reshape(greenFreqForRec, nLength, nRecs);
    %
    %         % calculate bigL matrix
    %         [meshXRec, meshXShot] = meshgrid(xRecGrid, xShotGrid);
    %         bigL = w(iw)^2 * (-rw1dFreq(iw)) * (greenFreqForShot(:, meshXShot(:)).' .* greenFreqForRec(:, meshXRec(:)).');
    %         % the true Hessian matrix 1
    %         hessianTrue = hessianTrue + bigL' * bigL;
    %         % the true Hessian matrix 2
    %         hessianTrue2 = hessianTrue2 + w(iw)^4 * abs(rw1dFreq(iw))^2 ...
    %             * ((conj(greenFreqForShot) * greenFreqForShot.') .* (conj(greenFreqForRec) * greenFreqForRec.'));
    %         mig = mig + bigL' * reshape(dataDeltaFreq(:, :, iw).', nShots * nRecs, 1);
    %         mig2 = mig2 + w(iw)^2 * (-conj(rw1dFreq(iw))) * sum(conj(greenFreqForShot) .* (conj(greenFreqForRec) * dataDeltaFreq(:, :, iw)), 2);
    %         % dm = (real(hessianTrue) \ real(rtm));
    %
    %         % u1 = w(iw)^2 * (-rw1dFreq(iw)) * (repmat(dm, 1, nShots) .* greenFreqForShot).' * greenFreqForRec;
    %         % u1 = bigL * dm;
    %         % u1 = reshape(u1, nShots, nRecs);
    %
    %         % for ixs = 1:nShots
    %         %     for ixr = 1:nRecs
    %         %         L = w(iw)^2 * (-rw1dFreq(iw)) * (greenFreqForShot(:, ixs).' .* greenFreqForRec(:, ixr).');
    %         %         u1_bak(ixs, ixr) = L * dm;
    %         %     end
    %         % end
    %
    %         % bias = u1 - dataDeltaFreq(:, :, iw).';
    %         % bias = reshape(bias, nShots * nRecs, 1);
    %         % gradVerify = gradVerify + real(bigL' * bias);
    %         % bias = reshape(bias, nShots, nRecs);
    %
    %         % for ixs = 1:nShots
    %         %     for ixr = 1:nRecs
    %         %         L = w(iw)^2 * (-rw1dFreq(iw)) * (greenFreqForShot(:, ixs).' .* greenFreqForRec(:, ixr).');
    %         %         gradVerify_bak = gradVerify_bak + ...
    %         %             real(L' * bias(ixs, ixr));
    %         %     end
    %         % end
    %
    %         timePerFreq = toc;
    %         fprintf('elapsed time = %fs\n', timePerFreq);
    %     end
    %
    %     lambda = max(diag(hessianTrue));
    %     hessianTrue = hessianTrue + lambda * eye(nLength, nLength);
    %
    %     % updated model
    %     modelOld = reshape(modelOld, nLength, 1);
    %     modelNew = modelOld + EPSILON * (real(hessianTrue) \ real(mig));
    %     modelOld = reshape(modelOld, nz, nx);
    %     modelNew = reshape(modelNew, nz, nx);
    %     modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    %     modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    %     vmNew = 1./sqrt(modelNew);
    %     vmNew = extBoundary(vmNew, nBoundary, 2);
    %
    %     figure(2);
    %     imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    %     xlabel('Distance (m)'); ylabel('Depth (m)');
    %     title('Updated Velocity Model');
    %     colormap(seismic);
    %     % save current updated velocity model
    %     % filenameVmNew = [pathVelocityModel, sprintf('/vmNew%d.mat', iter)];
    %     % save(filenameVmNew, 'vmNew', 'modelNew', '-v7.3');
    %
    %     dataDeltaFreq = zeros(nRecs, nShots, nfft);
    %% debug end
    
    clear('dataDeltaFreq');
    % load received surface data
    load(filenameDataTrueFreq);
    
    dataDeltaFreq = zeros(nRecs, nShots, nfft);
    
    % update dataDeltaFreq
    for ixs = 1:nShots %21:nx+20 % shot loop
        
        curXsPos = xShotGrid(ixs) + nBoundary; % shot position on x
        
        % generate shot signal
        sourceTime = zeros([size(V), nt]);
        sourceTime(1, curXsPos, :) = reshape(rw1dTime, 1, 1, nt);
        
        % generate shot record
        tic;
        [dataSmooth, ~] = fwdTimeCpmlFor2dAw(vmNew, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
        timeForward = toc;
        fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, elapsed time = %fs\n', curXsPos-nBoundary, x(curXsPos-nBoundary), timeForward);
        
        dataSmooth = dataSmooth(xr, :);
        
        dataDeltaFreq(:, ixs, :) = squeeze(dataTrueFreq(:, ixs, :)) - fftshift(fft(dataSmooth, nfft, 2), 2);
        
        % % debug begin
        % snapshotSmooth "blows out" when the wave propagates to the "abnormalities" in vmNew
        % figure(3);
        % for iframe = 1:nt
        %     imagesc(snapshotSmooth(1:end-nBoundary, nBoundary+1:end-nBoundary, iframe));
        %     title(['Iteration: ',num2str(iframe)])
        %     colorbar;
        %     drawnow
        % end
        % % debug end
        
    end % end shot loop
    filenameDataDeltaFreq = [pathVelocityModel, sprintf('/dataDeltaFreq%d.mat', iter)];
    save(filenameDataDeltaFreq, 'dataDeltaFreq', '-v7.3');
    
    % clear variables and functions from memory
    clear('dataTrueFreq');
    clear('dataDeltaFreq');
    clear('sourceTime');
    
    fprintf('Full-wave inversion iteration no. %d, model norm difference = %.6f\n', ...
        iter, norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro'));
    
    iter = iter + 1;
    
end

%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));
