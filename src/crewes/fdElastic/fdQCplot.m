function fdQCplot(modelDir)
    [parmFile] = fdFindParms(modelDir);
    %[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,....
    [Dt,Dxz,xMin,lengthX,lengthZ,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,....
        gfdF,wncvar,mvTif,contBlk] = readParmsB(parmFile);
        nx = round(lengthX/Dxz)+1;
        nz = round(lengthZ/Dxz)+1;
        %gfdFile] = readParmsB(parmFile);
        gfdFile = [modelDir,'\',gfdF];
        [ix1,iz1,nxf,nzf,nxplot,initzp,nzplot] =.....
    fdInitB3sizes(Dxz,nx,nz,lengthX,0,lengthZ,wncvar); %,iLbnd,shotX);
    %fdInitB3sizes(Dxz,nx,nz,mvXmax,mvZtop,mvZmax,wncMat); %,iLbnd,shotX);
        [Lp2m,Mu,Rho,LKrat,gfdFile,zMin] =...
    fdInitModB7(modelDir,gfdFile,Dxz,ix1,iz1,nxf,nzf,xMin,contBlk);
    Htitle = ['  Vp   ';'  Vs   ';'Density'];
    Vel(:,:,1) = sqrt(Lp2m./Rho);
    Vel(:,:,2) = sqrt(Mu./Rho);
    Vel(:,:,3) = Rho;
    vectX = [1:nxf]*Dxz;
    vectZ = [1:nzf]*Dxz;
    for iP = 1:3
        figure
        contourf(vectX,vectZ,Vel(:,:,iP)')
        title(Htitle(iP,:))
        axis ij
        colormap jet
        colorbar
        xlabel('X co-ordinate')
        ylabel('Depth')
        %boldlines
        bigfont(gca,1.5,2,1)
        whitefig
        grid off
    end
