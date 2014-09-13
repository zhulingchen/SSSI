function T = ray2d(V,Shot,dx)
%2D ray-tracing
%   RAY2D
load dA

%% Constants
v0 = 10000;
sz = 7;
sx = 7;
sz2 = 2*sz-1;
sx2 = 2*sx-1;
zs = Shot(1) + sz - 1;
xs = Shot(2) + sx - 1;

%% Derived values
V = 1./V;           % convert to slowness
mxV = max(V(:));
[nz,nx] = size(V);

%% Preallocate
T = ones(nz+sz2,nx+sx2) * v0;
M = T;
S = ones(nz+sz2-1,nx+sx2-1);

%%
M(sz:nz+sz,sx:nx+sx) = 0;
%dA=subs2(sz,sx);

iz = sz:nz+sz-1;
ix = sx:nx+sx-1;
S(iz,ix) = V;
S(nz+sz,ix)=2*S(nz+sz-1,ix)-S(nz+sz-2,ix);
S(iz,nx+sx)=2*S(iz,nx+sx-1)-S(iz,nx+sx-2);
S(nz+sz,nx+sx)=2*S(nz+sz-1,nx+sx-1)-S(nz+sz-2,nx+sx-2);

T(zs,xs) = 0;
M(zs,xs) = v0;

z1 = -sz+1:sz-1; z2 = -sz+1:sz-2; z3 = z1+zs;
x1 = -sx+1:sx-1; x2 = -sx+1:sx-2; x3 = x1+xs;
AS = S(z2+zs,x2+xs);
TT = T(z3,x3);
T(z3,x3) = min(reshape(dA*AS(:)+T(zs,xs),sz2,sx2),TT);
mxT = max(max(T(zs-1:zs+1,xs-1:xs+1)));

while true
    indx = T+M <= mxT+mxV;
    
    if isempty(indx)%idz)
         indx = M == 0;%[idz,idx] = find(M == 0);%break;
    end
    [idz,idx] = find(indx);
    M(indx) = v0;
    for i = 1:length(idz)
        z = idz(i);
        x = idx(i);
        mxT = max(mxT,T(z,x));
        AS = S(z+z2,x+x2);
        z3 = z+z1;
        x3 = x+x1;
        TT = T(z3,x3);
        T(z3,x3) = min(reshape(dA*AS(:)+T(z,x),sz2,sx2),TT);
    end %for
    
%     figure(3)
%     subplot(1,2,1)
%     imagesc(T*dx)
%     caxis([0 80])
%     subplot(1,2,2)
%     imagesc(M)
%     drawnow

    if prod(double(all(M(iz,ix))))
        break;
    end
    mxT = max(max(T(idz,idx)));
end %while

T = T(iz,ix)*dx;