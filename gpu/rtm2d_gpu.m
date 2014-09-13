function [model snapshot] = rtm2d_gpu(v,data,nz,dz,nx,dx,nt,dt)
%
% data(nx,nt)       data vector
% v(nz,nx)          velocity model
% nx                number of horizontal samples
% nz                number of depth samples
% nt                numer of time samples
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample

% add grid points for boundary condition

%v = [repmat(v(:,1),1,20), v, repmat(v(:,end),1,20)];
%v(end+20,:) = v(end,:);
%% Initialize storage
[nz,nx] = size(v);
[~,nt] = size(data);

%% Initialize GPU kernel and storage
k = parallel.gpu.CUDAKernel('rtm2d_kernel.ptx','rtm2d_kernel.cu');
k.ThreadBlockSize = [16,16,1];
fdm1 = gpuArray([data(:,nt)'; zeros(nz-1,nx)]);
fdm2 = gpuArray([data(:,nt-1)'; zeros(nz-1,nx)]);
fdm3 = gpuArray([data(:,nt-2)'; zeros(nz-1,nx)]);
data = gpuArray(data);

%% Boundary Absorbing Model
iz = 1:20;
boundary = gpuArray((exp(-( (0.015*(20-iz)).^2 ) )).^10);

%% Reverse-Time Migration
% finite difference coefficients
a = gpuArray((v*dt/dx).^2);    % wave equation coefficient
b = 2-4*a;

cz = 3;
snapshot = zeros(nz,nx,nt);
for it = (nt-1):-1:1
    cz = cz+1;
    bz = min(cz,nz);
     
    % call gpu function
    [fdm1,fdm2,fdm3] = feval(k,fdm1,fdm2,fdm3,boundary,...
        a,b,nz,nx,bz,it,data,nt);

    %
    tmp = gather(fdm1);
    imagesc(tmp)
    title(['Time index: ',num2str(it)])
    drawnow
    %}
    snapshot(:,:,it) = tmp;
end % time loop

% write out final wavefield
model = gather(fdm1);
