function [wncmat,iZlvl,thrat].....
    = fdInitB2Gen(mvPlot,Dt,Dxz,wncvar)
%function [wncmat,thrat].....
%    = fdInitB2Gen(mvPlot,Dt,Dxz,wncvar)
% function [nxplot,initzp,nzplot,wncmat,thrat].....
%     = fdInitB2Gen(mvXmax,mvZtop,mvZmax,mvPlot,Dt,Dxz,wncvar)
%Initialize general conditions for an FD model
%
%The input parameters are
% %mvXmax  .... X-length of movie frames to display or plot
% %mvZtop  .... Z-level at top of each movie frame, often 0.
% %mvZmax  .... Z-depth of movie frames to display or plot
%mvPlot  .... The code number of the movie snapshot plot
%Dt      .... FD sample rate in seconds
%Dxz     .... FD sample rate in (metres)
%wncvar  .... The wavenumber correction file (within quotes), ('' is none) 
%
%The output parameters are
% %nxplot  .... No. of X points to plot
% %initzp  .... Initial Z point to plot
% %nzplot  .... No. of Z points to plot
%wncmat  .... A set of FD correction matrices, for particular Vp and Vs
%thrat   .... (Dt/Dxz)^2, used in the FD calculations

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

% nxplot = round(mvXmax/Dxz)+1; 
% initzp = round(mvZtop/Dxz)+1; 
% nzplot = round(mvZmax/Dxz)+1;
%disp([nxplot,initzp,nzplot])
%wnc = [];
%nShot = 1;      %For now
% wnc = 1;
% icorr = 0;
% load wncUnity
% wncU = wnc;
% %disp(wnc)
% load wncEdge
% %disp(wnc)
%supp = 0;
% if iTbnd > 0
%     supp = suppress;                %Set up supp variable here
% end
wncmat = 1; iZlvl = 1;
%small = 0.0000001;   %0.00001;
if ~isempty(wncvar)
    loadStr = ['load ' wncvar];
    disp(loadStr)
    eval(loadStr)       %Corrections in variable 'wncmat'
    disp(size(wncmat)); disp(size(iZlvl));
    
    %disp(wncmat)
end
% figure
% flipy
set(gca,'NextPlot','replacechildren')
%jfr = 0;
% %mint = round((nstep-mvTinit)/nframe);
% mint = round(nstep/nframes);
if mvPlot>=0
    set(gca,'xlimmode','manual')
end
thrat = (Dt/Dxz)^2;

