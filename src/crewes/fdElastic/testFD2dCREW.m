%testFD2d - Run 2D finite-difference model
clear
global tictoc
tictoc = clock;
 tic
 
% F2dmats6                           %Latest
% F2dmats7                           %Latest

% load('c3p20s10m50.mat')
% whos

 surfB5FD2d    %compB4AllFD2d     Getr4F resampSFplot %Latest

% [generic] = fdReadgeo2('zone.geo');

% [Dt,Dxz,fractx,nWide,iPlTerm,wncvar,pvel,svel] = ....
%            readParmsCr4('test.parc');

%disp(Uz(1:3,1:3))
%disp(Ux(1:3,1:3))
% figure
% plot(surfUz(2,:))
% grid on
% title('Corrected, ixm2 = 2:ix9')
% print -dtiff -r150 corr2.tif
% xA = 1:400;
% plot(xA,wellUx(140,:),xA,wellUx(30,:))
% grid on
%help deconw
%help seismic
%help ormsby
%help resamp
%help sinci
%help plotseis
%help displaytools
%help plotimage
% figure
% whitefig
% print -dtiff -r150 whiteF.tif
%whos
%help fftrl
 toc
