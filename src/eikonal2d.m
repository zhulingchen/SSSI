function T = eikonal2d(V, dx, sz, sx, iMax)
% Solution of eikonal equation in 2D rectangular domain for a single source
% at (sz,sx).
%
% V                 velocity model
% dx                grid spacing, assume dx=dy
% sx                x coordinates of the seismic events
% sz                z coordinates of the seismic events 
% This matlab source file is free for use in academic research.
% iMax              maximum number of iterations of sweeping (related to the size of model and initial condition)
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology
%
%Reference: H. Zhao, A fast sweeping method for Eikonal equations,
%Mathematics of computation, 74(250), pp. 603-627, 2004

%% initialization 
[nz, nx]= size(V);
Slow = 1./V;
Told = 1e6*ones(nz,nx);
Tnew = 1e6*ones(nz,nx);
Told(sz,sx)=0;


for iter=1:iMax
%% first Sweep i=1:nz, j=1:nx
for i=1:nz
    for j=1:nx
      if (i==sz && j==sx)
         Tnew(i,j) = Told(i,j);
      else
          
         if i==1
             a=Told(2,j);
         elseif i==nz
             a=Told(nz-1,j);
         else
             a=min(Told(i-1,j),Told(i+1,j));
         end
         
         if j==1
             b=Told(i,2);
         elseif j==nx
             b=Told(i,nx-1);
         else
             b=min(Told(i,j-1),Told(i,j+1));
         end
         
         if (abs(a-b)<Slow(i,j)*dx)
             TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
         else
             TT=min(a,b)+Slow(i,j)*dx;
         end
         Tnew(i,j) = min(Told(i,j),TT);
      end
   end
end

Told=Tnew;

%% second sweep
for i=nz:-1:1
    for j=1:nx
      if (i==sz && j==sx)
         Tnew(i,j) = Told(i,j);
      else
          
         if i==1
             a=Told(2,j);
         elseif i==nz
             a=Told(nz-1,j);
         else
             a=min(Told(i-1,j),Told(i+1,j));
         end
         
         if j==1
             b=Told(i,2);
         elseif j==nx
             b=Told(i,nx-1);
         else
             b=min(Told(i,j-1),Told(i,j+1));
         end
         
         if (abs(a-b)<Slow(i,j)*dx)
             TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
         else
             TT=min(a,b)+Slow(i,j)*dx;
         end
         Tnew(i,j) = min(Told(i,j),TT);
      end
   end
end

Told=Tnew;

%% Third sweep
for i=nz:-1:1
    for j=nx:-1:1
      if (i==sz && j==sx)
         Tnew(i,j) = Told(i,j);
      else
          
         if i==1
             a=Told(2,j);
         elseif i==nz
             a=Told(nz-1,j);
         else
             a=min(Told(i-1,j),Told(i+1,j));
         end
         
         if j==1
             b=Told(i,2);
         elseif j==nx
             b=Told(i,nx-1);
         else
             b=min(Told(i,j-1),Told(i,j+1));
         end
         
         if (abs(a-b)<Slow(i,j)*dx)
             TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
         else
             TT=min(a,b)+Slow(i,j)*dx;
         end
         Tnew(i,j) = min(Told(i,j),TT);
      end
   end
end

Told=Tnew;

%% Fourth sweep
for i=1:nz
    for j=nx:-1:1
      if (i==sz && j==sx)
         Tnew(i,j) = Told(i,j);
      else
          
         if i==1
             a=Told(2,j);
         elseif i==nz
             a=Told(nz-1,j);
         else
             a=min(Told(i-1,j),Told(i+1,j));
         end
         
         if j==1
             b=Told(i,2);
         elseif j==nx
             b=Told(i,nx-1);
         else
             b=min(Told(i,j-1),Told(i,j+1));
         end
         
         if (abs(a-b)<Slow(i,j)*dx)
             TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
         else
             TT=min(a,b)+Slow(i,j)*dx;
         end
         Tnew(i,j) = min(Told(i,j),TT);
      end
   end
end

Told=Tnew;
end


T=Told;

