function [Lp2m,Mu,Rho,LKrat,oldFile,zMin] =...
    fdInitModB6(oldFile,Dxz,ix1,iz1,nxf,nzf,xMin,contBlk)
% function [Lp2m,Mu,Rho,LKrat,oldFile,zMin] =...
%     fdInitModB6(oldFile,Dxz,nx,nz,ix1,iz1,nxf,nzf,xMin,contBlk)
%Get the geological model from disc - fill in FD grid
%The input parameters are
%oldFile .... The geological definition file (within quotes) (gfdfile)
%Dxz     .... FD sample rate in (metres)
%nx
%nz
%ix1   ...... X co-ordinate of start of model in arrays
%iz1   ...... Z co-ordinate of start of model in arrays
%nxf     .... Number of X samples including border
%nzf     .... Number of Z samples including border
%xMin    .... Left edge of the FD model extracted from the geological model
%contBlk .... The indicator of a 'cont' (continuous) or 'block' (blocked)
%               geological model
%The output parameters are
%Lp2m    .... Lamda plus 2*Mu
%Mu      .... Mu
%Rho     .... Density
%LKrat   .... (Lp2m-2*Mu)./Lp2m
%oldFile .... The geological definition file (within quotes) (gfdfile)
%zMin    .... Top of geological model
%
% P.M. Manning, Dec 2011
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

if(~strncmp(contBlk,'log',3))       %If 'log' go to bottom of listing
    UpDown = 'Down';
    %contBlk = 'cont'; %''; %'cont';
%     if isempty(oldFile)
%         disp('The geology definition files for F-D input are ')
%         ls *.geo
%         [oldFile] = fdLoadModelSt;  %??
%     end
    disp(oldFile)
    [generic] = fdReadgeo2(oldFile);
    load(generic,'-mat')
%     whos
%     disp(Vp);disp(Vs);disp(rho);
%     disp(Xtops);disp(Ztops);disp(nEnt);
%     stop
Lp2m = zeros(nxf,nzf);
Mu = zeros(nxf,nzf);
Rho = zeros(nxf,nzf);
Lp2mH = rho.*Vp.^2;     %Take dimensions of input vectors
MuH = rho.*Vs.^2;

%nxu = nxf-1;
%nzu = nzf-1;
nxu = nxf;
nzu = nzf;

%Fill whole area with top values
Lp2m(1:nxf,1:nzf) = ones(nxf,nzf)*Lp2mH(1);
Mu(1:nxf,1:nzf) = ones(nxf,nzf)*MuH(1);
Rho(1:nxf,1:nzf) = ones(nxf,nzf)*rho(1);

Wons = ones(1,nzu);
nTops = size(Xtops);
disp(Vp)
%nTops = size(Vp);
zMin = Ztops(1,1);
%whos
Lp2mG = zeros(size(Vp));
MuG = Lp2mG; RhoG = Lp2mG;
%disp(size(Lp2m)); disp(size(Ztops)); 
if strncmpi('cont',contBlk,4)
    for iTz = 1:nTops-1
        %Lp2mG(iTz) = (Vp(iTz+1).^2-Vp(iTz).^2)*rho(1,iTz)/(Ztops(1,iTz+1)-Ztops(1,iTz));
        
        %Lp2mG(iTz) = (Vp(iTz+1).^2-Vp(iTz).^2)*rho(iTz,1)/(Ztops(iTz+1,1)-Ztops(iTz,1));
        %MuG(iTz) = (Vs(iTz+1).^2-Vs(iTz).^2)*rho(iTz,1)/(Ztops(iTz+1,1)-Ztops(iTz,1));
        %RhoG(iTz) = (rho(iTz+1)-rho(iTz))/(Ztops(iTz+1,1)-Ztops(iTz,1));
        %Zdist = (Ztops(iTz+1,1)-Ztops(iTz,1))/Dxz;
        Zdist = Ztops(iTz+1,1)-Ztops(iTz,1);
        %disp(Zdist)
%         Lp2mG(iTz) = (Vp(iTz+1).^2-Vp(iTz).^2)*rho(iTz,1)/Zdist;
%         MuG(iTz) = (Vs(iTz+1).^2-Vs(iTz).^2)*rho(iTz,1)/Zdist;
        Lp2mG(iTz) = (Vp(iTz+1).^2*rho(iTz+1)-Vp(iTz).^2*rho(iTz))/Zdist;
        MuG(iTz)   = (Vs(iTz+1).^2*rho(iTz+1)-Vs(iTz).^2*rho(iTz))/Zdist;
        %MuG(iTz) = (Vs(iTz+1).^2-Vs(iTz).^2)*rho(iTz,1)/Zdist;
        RhoG(iTz) = (rho(iTz+1)-rho(iTz))/Zdist;
    end
end
%disp(Lp2mG)
%for iT=2:nTops              %Go from the second row
for iT=1:nTops              %Go from the first row (no elevation profile)
    xT = xMin;
    for ix=ix1:nxu
        ig = find(Xtops(iT,:)>xT);
        zTp = Ztops(iT,ig(1))+(xT-Xtops(iT,ig(1)))*(Ztops(iT,ig(1))-Ztops(iT,ig(1)-1)).....
            /(Xtops(iT,ig(1))-Xtops(iT,ig(1)-1));
        iTz = ceil(zTp/Dxz)+iz1;    %Index of first sample below top iT
                                    %New velocities from here
        ziTz = (iTz-iz1)*Dxz;
        %disp(iTz)
        %stop
%         if strcmpi(UpDown,'Up')
%             if iTz>1
%                 if iTz>nz; iTz = nz; end
%                 Lp2m(ix,1:iTz) = Lp2mH(iT)*Wons(1:iTz);
%                 Mu(ix,1:iTz) = MuH(iT)*Wons(1:iTz);
%                 Rho(ix,1:iTz) = rho(iT)*Wons(1:iTz);
%             end
%         else
            %iTz = iTz+1;
            %if iTz<nz
            if iTz<nzf
                %if iTz<1; iTz = 1; end
                if ~strcmpi('cont',contBlk)         %Blocked
                    Lp2m(ix,iTz:nzu) = Lp2mH(iT)*Wons(iTz:nzu);
                    Mu(ix,iTz:nzu) = MuH(iT)*Wons(iTz:nzu);
                    Rho(ix,iTz:nzu) = rho(iT)*Wons(iTz:nzu);
                else                                %Continuous
                    %izA = 0;
                    %disp(RhoG(1));disp(Dxz);disp(rho(1))
                    for iz = iTz:nzu
%                         Lp2m(ix,iz) = Lp2mH(iT)+izA*Lp2mG(iT);
%                         Mu(ix,iz) = MuH(iT)+izA*MuG(iT);
%                         Rho(ix,iz) = rho(iT)+izA*RhoG(iT);
%                         izA = izA+1;
                        Lp2m(ix,iz) = Lp2mH(iT)+(ziTz-zTp)*Lp2mG(iT);
                        Mu(ix,iz) = MuH(iT)+(ziTz-zTp)*MuG(iT);
                        Rho(ix,iz) = rho(iT)+(ziTz-zTp)*RhoG(iT);
                        ziTz = ziTz+Dxz;
%                         disp(Rho(1,iz))
%                         pause
                    end
                    %stop
                end
            end
%         end
        xT = xT+Dxz;
    end
    %disp(iTz)
end
%Fill in parameters at top
for iz = 1:iz1-1
    Lp2m(:,iz) = Lp2m(:,iz1+1);
    Mu(:,iz) = Mu(:,iz1+1);
    Rho(:,iz) = Rho(:,iz1+1);
end
zVect = 1:nzf;
iXpos = round(nxf*0.5);
%disp(Lp2m(iXpos,1:9))
figure
plot(zVect,Lp2m(iXpos,:),zVect,Mu(iXpos,:))
%disp('Lp2m'); disp(Lp2m(iXpos,nzf-9:nzf));
%disp('Mu'); disp(Mu(iXpos,nzf-9:nzf));

title('Lp2m and Mu profiles at centre')
%plot(zVect,Rho(iXpos,:))
%disp(Rho(1,1:9))
%stop
else                %This 'log' option models space from one log
    load(oldFile)
    nzLog = length(vpblk);
    %disp(vpblk(nzLog-9:nzLog))      %Likely highest velocities
    avP = sum(vpblk)/nzLog;
    avS = sum(vsblk)/nzLog;
    disp(['Average P = ',num2str(avP)])
    disp(['Average S = ',num2str(avS)])
    VpH = ones(1,nzf)*vpblk(nzLog);
    VsH = ones(1,nzf)*vsblk(nzLog);
    rhoH = ones(1,nzf)*rhoblk(nzLog);
    if nzLog>nzf
        nzLog = nzf;
    end
    rhoH(1:nzLog) = rhoblk(1:nzLog);
    VpH(1:nzLog) = vpblk(1:nzLog);
    VsH(1:nzLog) = vsblk(1:nzLog);
    Lp2mH(1:nzLog) = rhoblk(1:nzLog).*VpH.^2;
    MuH(1:nzLog) = rhoblk(1:nzLog).*VsH.^2;
    Wons = ones(nxf,1);
    Rho = Wons*rhoH;
    Lp2m = Wons*Lp2mH;
    Mu = Wons*MuH;
    zMin = zblk(1);
end
%figure
%plot(Lp2m(1,:))
LKrat = (Lp2m-2*Mu)./Lp2m;
%disp('HERE')
