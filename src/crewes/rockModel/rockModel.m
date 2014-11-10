clear;
home;
% define the uncracked rock's properties
Vp1=3.026;
Vs1=1.721;
Dens1=2.133;
% define parameters for inclusions
PHI=0.005;
alpha=0.01;
Dens2=1.1;
Vp2=1.43;
Vs2=0;


[Vp_KT,Vs_KT]=crack_KT(Dens1,Vp1,Vs1,PHI,alpha,Dens2,Vp2,Vs2,1); % spheroid = 1, sphere =0
[Vp_KTB,Vs_KTB]=crack_KTB(Dens1,Vp1,Vs1,PHI,alpha,Dens2,Vp2,Vs2,4); % penny crack = 4
[Vp_Hud,Vs_Hud]=crack_Hudson(Dens1,Vp1,Vs1,PHI,alpha,Dens2,Vp2,Vs2,1); % penny crack with fluid substitution = 1

