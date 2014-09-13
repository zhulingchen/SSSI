%%
clc, close all, clear all
load dat4gpu
dt = 0.0020;
!nvcc -ptx fm2d_kernel.cu
tic, [datag snapshotg] = fm2d_gpu(V,rw,nz,dz,nx,dx,nt,dt); toc
%%
close all, clc
!nvcc -ptx rtm2d_kernel.cu
tic, [modelg rtmsnapshotg] = rtm2d_gpu(V,data,nz,dz,nx,dx,nt,dt); toc;

%%
show = true;
[nz,nx,nt]=size(rtmsnapshotg);

M = 0;
s2 = 0;
for i = 1:nt
    M = snapshotg(:,:,i).*rtmsnapshotg(:,:,nt-i+1)+M;
    s2 = snapshotg(:,:,i).^2+s2;
    
    if show
        subplot(2,2,3)
        imagesc(snapshotg(1:end-20,21:end-20,i))
        
        subplot(2,2,4)
        imagesc(rtmsnapshotg(1:end-20,21:end-20,nt-i+1))
        
        subplot(2,2,2)
        imagesc(diff(M(1:end-20,21:end-20)./s2(1:end-20,21:end-20),2,1))
        
        drawnow
    end
end