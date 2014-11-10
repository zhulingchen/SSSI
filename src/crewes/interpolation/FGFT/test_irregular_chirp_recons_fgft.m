% Testing forward fgft
% Author: Mostafa Naghizadeh; Copyright (C) 2010
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

clc
%close all;
clear all;

% Chirp example
dt = 4/1000;,
N = 256;
% Compute a test signal (a hyperbolic chirp)
[d,t,cf] = mychirp(10,110,dt,N);
% Add noise to the chirp
%d = d + 0.6*randn(size(d));
%d=[d;flipud(d)];
N=256;

% load fx_slice.mat
% d=x1.';
d_orig=d;

% Irregular interpolation
% Randomizing the offset locations	
kn1=randperm(N);
% Randomly eliminating some traces
percent=0.5;
kn2=kn1(1:ceil(percent*N));
% sorting random chosen traces in ascending order
h=[sort(kn2)];
%h=[1:2:N];  % Regular decimation
%h=[90:100];  % Gap
%save h_rand.mat h;
%load h_rand.mat
d(h)=0;

%Finding the indices of known components
h1=[1:N];
[a,b]=wcommon(h,h1);
H=find(b-1);

% Forward FGFT
D=forward_fgft(d');

% Tresholding smaller values of D
D1=forward_fgft(d_orig');
pctg=0.8;
cfs = sort(abs(D1));
nb = round(pctg*length(cfs));
cutoff = cfs(nb);
% Set small coefficients to zero
for w=1:length(D1)
    MASK(w)= (abs(D1(w))>cutoff);
end

D=scale_fgft(D);
PK=plot_fgft(D);

d3=d(H);
M=MASK;
nh=N;
iter_cg=3;
iter_bl=50;

% For irregular sampling
INTD=irls_fitting_fgft(d3.',H,nh,iter_cg,iter_bl);

% % For regular sampling provided that we know the mask
% INTD=ls_mask_fitting_fgft(d3.',H,M,nh,iter_cg);

%INTD=INTD/max(abs(INTD));
figure
subplot(411)
h2=plot(linspace(1,N,N),d,'k');axis tight;
setit(h2,'Time samples','Amplitude');

subplot(412)
colormap(flipud(gray));
h2=imagesc([1:N],linspace(0,.51,129),abs(PK(:,:)));
setit02(h2,'Time samples','Frequency');

subplot(413)
D2=forward_fgft(INTD);
D2=scale_fgft(D2);
PK3=plot_fgft(D2);

colormap(flipud(gray));
h2=imagesc([1:N],linspace(0,.51,129),abs(PK3(:,:)));
setit02(h2,'Time samples','Frequency');

subplot(414)
h2=plot(linspace(1,N,N),real(INTD),'k'); axis tight;
setit(h2,'Time samples','Amplitude');
