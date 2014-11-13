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

assert(isa(V, 'double'));
assert(isa(dx, 'double'));
assert(isa(sz, 'double'));
assert(isa(sx, 'double'));
assert(isa(iMax, 'double'));

%% initialization
[nz, nx]= size(V);
Slow = 1./V;
Slow = padarray(Slow, [1, 1], 'replicate');
Told = 1e6 * ones(nz+2,nx+2);
Tnew = 1e6 * ones(nz+2,nx+2);
Told(sz+1,sx+1)=0;
T = zeros(nz+2, nx+2);

a = 0;
b = 0;

for iter=1:iMax
    %% first Sweep i=1:nz, j=1:nx
    for i=2:nz+1
        for j=2:nx+1
            a=min(Told(i-1,j),Told(i+1,j));
            
            b=min(Told(i,j-1),Told(i,j+1));
            
            if (abs(a-b)<Slow(i,j)*dx)
                TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
            else
                TT=min(a,b)+Slow(i,j)*dx;
            end
            Tnew(i,j) = min(Told(i,j),TT);
        end
    end
    
    Tnew(sz+1,sx+1)=0;
    Told=Tnew;
    
    %% second sweep
    for i=nz+1:-1:2
        for j=2:nx+1
            a=min(Told(i-1,j),Told(i+1,j));
            
            b=min(Told(i,j-1),Told(i,j+1));
            
            if (abs(a-b)<Slow(i,j)*dx)
                TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
            else
                TT=min(a,b)+Slow(i,j)*dx;
            end
            Tnew(i,j) = min(Told(i,j),TT);
            
        end
    end
    Tnew(sz+1,sx+1)=0;
    Told=Tnew;
    
    %% Third sweep
    for i=nz+1:-1:2
        for j=nx+1:-1:2
            
            a=min(Told(i-1,j),Told(i+1,j));
            
            b=min(Told(i,j-1),Told(i,j+1));
            
            if (abs(a-b)<Slow(i,j)*dx)
                TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
            else
                TT=min(a,b)+Slow(i,j)*dx;
            end
            Tnew(i,j) = min(Told(i,j),TT);
            
        end
    end
    Tnew(sz+1,sx+1)=0;
    Told=Tnew;
    
    %% Fourth sweep
    for i=2:nz+1
        for j=nx+1:-1:2
            
            a=min(Told(i-1,j),Told(i+1,j));
            
            b=min(Told(i,j-1),Told(i,j+1));
            
            if (abs(a-b)<Slow(i,j)*dx)
                TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
            else
                TT=min(a,b)+Slow(i,j)*dx;
            end
            Tnew(i,j) = min(Told(i,j),TT);
            
        end
    end
    Tnew(sz+1,sx+1)=0;
    Told=Tnew;
end


T=Told(2:end-1, 2:end-1);

