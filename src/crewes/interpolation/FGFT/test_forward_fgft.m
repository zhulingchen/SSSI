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
N =256;
% Compute a test signal (a hyperbolic chirp)
[d,t,cf] = mychirp(10,115,dt,N);
% Add noise to the chirp
%d = d + 0.4*randn(size(d));
%d=[d;flipud(d)];
N=256;

% % Box car example
% n=256;
% l=30;
% BX=zeros(1,n);
% BX(floor(n/2-l/2)+1:floor(n/2+l/2))=1;
% d=BX';

% % Sin example
% x1 = linspace(-50*pi,50*pi,256);
% x2=linspace(-10*pi,10*pi,256);
% d1=sin(x1); d2=sin(x2); d=d1+d2;
% d(1:10)=d(1:10).*(1./fliplr(linspace(1,4,10)));
% d(247:256)=d(247:256).*(1./linspace(1,4,10));
% d=d';

% % Irregular interpolation
% % Randomizing the offset locations	
% kn1=randperm(N);
% % Randomly eliminating some traces
% percent=0.3;
% kn2=kn1(1:ceil(percent*N));
% % sorting randomly chosen traces in ascending order
% h=[sort(kn2)];
% h=[120:170];
% d(h)=0;

% % Missing parts of d
% % Sampling operator is T 
%  T = ones(size(d));
%  T(20:30) = 0;
%  T(50:60) = 0;
%  T(120:147) = 0;
%  T(190:197) = 0;
% 
% % Sample the data
%  d = d.*T;

D=forward_fgft(d');
% D=smooth(D,5).';

% % Tresholding smaller values of D
% pctg=0.75;
% cfs = sort(abs(D));
% nb = round(pctg*length(cfs));
% cutoff = cfs(nb);
% % Set small coefficients to zero
% for w=1:length(D)
%     if (abs(D(w))<cutoff)
%         D(w)=0;
%     end
% end

D1=scale_fgft(D);

PK=plot_fgft(D1);


figure
subplot(411)
h2=plot(linspace(1,N,N),d,'k');axis ([1 256 -1.1 1.1]);
setit(h2,'Time samples','Amplitude');


subplot(412)
h2=plot(linspace(1,128,128),abs(D1(129:256))/max(abs(D1(129:256))),'k');axis tight;
setit(h2,'FGFT samples','Amplitude');


subplot(413)
colormap(flipud(gray));
h2=imagesc([1:N],linspace(0,.5,129),abs(PK(1:129,:)));
setit02(h2,'Time samples','Frequency');

X2=inverse_fgft(D);

subplot(414)
h2=plot(linspace(1,N,N),real(X2),'k'); axis ([1 256 -1.1 1.1]);
setit(h2,'Time samples','Amplitude');
