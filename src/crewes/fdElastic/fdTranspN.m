function [D2UxDx,D2UzDx,D2UxDz,D2UzDz,D2UxDxz,D2UzDxz] = fdTranspN(Ux,Uz,Lp2m,Mu,LKrat)
%function [D2UxDx,D2UzDx,D2UxDz,D2UzDz,D2UxDxz,D2UzDxz] = fdTransp(Ux,Uz,Lp2m,Mu,LKrat)
%Calculates the first usable layer of finite-difference second derivatives

[nxf,nzf] = size(Ux);
ix9 = nxf-1;
DUxDx = zeros(nxf,2);
DUxDz = DUxDx; DUzDx = DUxDx; DUzDz = DUxDx;
D2UxDx = zeros(nxf,1);
D2UzDx = D2UxDx; D2UxDxz = D2UxDx; D2UzDxz = D2UxDx; 
for iz = 1:nzf                % Loop through all Z values
  %Calculate derivatives of Ux with respect to x  (1)
    for ix = 2:nxf
        DUxDx(ix,iz) = (Ux(ix,iz)-Ux(ix-1,iz)).*Lp2m(ix,iz);
    end
  %Calculate derivatives of Uz with respect to x  (8)
    for ix = 1:ix9
        DUzDx(ix,iz) = (Uz(ix+1,iz)-Uz(ix,iz)).*Mu(ix,iz);
    end
end
iz = 2;
for ix = 2:ix9
    D2UxDx(ix) = DUxDx(ix+1,iz)-DUxDx(ix,iz);
end
for ix = 2:ix9
    D2UzDx(ix) = DUzDx(ix,iz)-DUzDx(ix-1,iz);
end
for iz = 2:nzf
    %Calculate derivatives of Ux with respect to z  (3)
    DUxDz(:,iz) = (Ux(:,iz)-Ux(:,iz-1)).*Mu(:,iz);
    %Calculate derivatives of Uz with respect to z  (7)
    DUzDz(:,iz-1) = (Uz(:,iz)-Uz(:,iz-1)).*Lp2m(:,iz); %1 shift
end
iz = 2;
for ix = 1:ix9
    D2UxDxz(ix+1) = (DUxDz(ix+1,iz)-DUxDz(ix,iz))...
    +(DUxDx(ix+1,iz)-DUxDx(ix+1,iz-1)).*LKrat(ix+1,iz);
    D2UzDxz(ix) = (DUzDx(ix,iz+1)-DUzDx(ix,iz))...
    +(DUzDz(ix+1,iz)-DUzDz(ix,iz)).*LKrat(ix,iz);
end
iz = 2;
D2UxDz = DUxDz(:,iz+1)-DUxDz(:,iz);
D2UzDz = DUzDz(:,iz)-DUzDz(:,iz-1);
%disp(D2UxDx(130:135))
% end

