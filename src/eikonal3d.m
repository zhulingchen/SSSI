function T = eikonal3d(V, dx, sz, sx, sy, iMax)
% Solution of eikonal equation in 2D rectangular domain for a single source
% at (sz,sx).
%
% V                 velocity model
% dx                grid spacing, assume dx=dy=dz
% sx                x coordinates of the seismic events
% sy                x coordinates of the seismic events 
% sz                z coordinates of the seismic events 
% iMax              maximum number of iterations of sweeping (related to the size of model and initial condition)
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology
%
%Reference: H. Zhao, A fast sweeping method for Eikonal equations,
%Mathematics of computation, 74(250), pp. 603-627, 2004

%% initialization 
[nz, nx, ny]= size(V);
Slow = 1./V; 
Told = 1e6*ones(nz,nx,ny);
Tnew = 1e6*ones(nz,nx,ny);
Told(sz,sx,sy) = 0; % set the source point


for iter=1:iMax
%% first Sweep i=1:nz, j=1:nx, k=1:ny
for i=1:nz
    for j=1:nx
        for k=1:ny
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end

Tnew(sz,sx,sy)=0;
Told=Tnew;

%% second Sweep i=nz:1, j=1:nx, k=1:ny
for i=nz:-1:1
    for j=1:nx
        for k=1:ny
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told=Tnew;


%% third Sweep i=nz:1, j=nx:1, k=1:ny
for i=nz:-1:1
    for j=nx:-1:1
        for k=1:ny
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told=Tnew;

%% fourth Sweep i=nz:1, j=nx:1, k=ny:1
for i=nz:-1:1
    for j=nx:-1:1
        for k=ny:-1:1
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told=Tnew;

%% fifth Sweep i=1:nz, j=nx:1, k=ny:1
for i=1:nz
    for j=nx:-1:1
        for k=ny:-1:1
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
            
       end
   end
end
Tnew(sz,sx,sy)=0;
Told=Tnew;

%% sixth Sweep i=1:nz, j=1:nx, k=ny:1

for i=1:nz
    for j=1:nx
        for k=ny:-1:1
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told=Tnew;

%% seventh Sweep i=1:nz, j=nx:1, k=1:ny

for i=1:nz
    for j=nx:-1:1
        for k=1:ny
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told = Tnew;

%% eighth Sweep i = nz:1, j=1:nx, k= ny:1

for i=nz:-1:1
    for j=1:nx
        for k=ny:-1:1
            
                
              if i==1
                   a=Told(2,j,k);
              elseif i==nz
                   a=Told(nz-1,j,k);
              else
                   a=min(Told(i-1,j,k),Told(i+1,j,k));
              end
         
                 if j==1
                     b=Told(i,2,k);
                 elseif j==nx
                     b=Told(i,nx-1,k);
                 else
                     b=min(Told(i,j-1,k),Told(i,j+1,k));
                 end
                 
                 if k==1
                     c=Told(i,j,2);
                 elseif k==ny
                     c=Told(i,j,ny-1);
                 else
                     c=min(Told(i,j,k-1),Told(i,j,k+1));
                 end
                 
                  TT = sqSolver(a,b,c,Slow(i,j,k)*dx);
                  Tnew(i,j,k) = min(Told(i,j,k),TT);
             
       end
   end
end
Tnew(sz,sx,sy)=0;
Told = Tnew;
end
T=Told;



