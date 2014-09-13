function [data snapshot] = fm2d_gpu(v,model,nz,dz,nx,dx,nt,dt)
%
% model(nz,nx)      model vector
% v(nz,nx)          velocity model
% nx                number of horizontal samples
% nz                number of depth samples
% nt                numer of time samples
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample

% add grid points for boundary condition
%model = [repmat(model(:,1),1,20), model, repmat(model(:,end),1,20)];
%model(end+20,:) = model(end,:);

%v = [repmat(v(:,1),1,20), v, repmat(v(:,end),1,20)];
%v(end+20,:) = v(end,:);
%% Initialize storage
[nz,nx] = size(model);
data = zeros(nx,nt);
snapshot = zeros(nz,nx,nt);

%% Initialize GPU kernel and storage
k = parallel.gpu.CUDAKernel('fm2d_kernel.ptx','fm2d_kernel.cu');
k.ThreadBlockSize = [16,16,1];
fdm1 = gpuArray(zeros(nz,nx));
fdm2 = gpuArray(model);
fdm3 = gpuArray(zeros(nz,nx));

%% Boundary Absorbing Model
iz = 1:20;
boundary = gpuArray((exp(-( (0.015*(20-iz)).^2 ) )).^10);

%% Forward-Time Modeling
data(:,1)  = model(1,:);

% finite difference coefficients
a = gpuArray((v*dt/dx).^2);
b = 2-4*a;

for it = 2:nt
    
    % call gpu function
    [fdm1,fdm2,fdm3] = feval(k,fdm1,fdm2,fdm3,boundary,a,b,nz,nx);
    % update data
    tmp = gather(fdm3);
    data(:,it) = tmp(1,:);
    
    %
    imagesc(tmp)
    title(['Time index: ',num2str(it)])
    drawnow
    
    snapshot(:,:,it) = tmp;
end % time loop

%data = data(21:nx-20,:);