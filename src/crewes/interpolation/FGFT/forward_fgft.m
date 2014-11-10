function K=forward_fgft(X)
% This function computes the Generalized Fourier transform based on a
% paper by Robert Brown et al. (2010)
% [IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 58, NO. 1, JANUARY 2010]
% The idea is to perform time-frequency analysis (Gabor, S-Transform,..)
% in a very efficient way, and also saving fewer coefficients (namely the
% same size as the original data).
%
%% Input
%   X: 1D signal (length of X should be a power of 2)
%
%% Output
%   Y: 1D general Fourier transform
%
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


% Determining the size of data
nx=length(X);

% Going to the next power of 2
n2=nextpow2(nx);

% Taking Fourier transform of input
K=fft(X);

% % First and last 2 samples
% WW=zeros(1,nx);
% WW=W;
% FWW=fftshift(fft(WW));
% % Choosing only the middle part of the Fourier window
% WFWW=FWW(nx/2:nx/2+1);
% WFWW=WFWW/max(abs(WFWW));
% 
% % Positive frequencies
% K(nx-1:nx)=(1/sqrt(2^(2-1)))*fft((K(nx-1:nx).*abs(WFWW)));
% K(1:2)=(1/sqrt(2^(2-1)))*fft((K(1:2).*abs(WFWW)));

% Starting loop to pick parts of K and apply FFT in smaller sizes
for ia=2:n2-1
    % At each ia step we will take the fft of positive and
    % negative frequencies with variant window size.
    pfs=2^ia-2^(ia-1)+1;
    pfe=2^ia;
    nfs=(2^n2-(2^ia-2^(ia-1)+1)+1)-2^(ia-1)+1;
    nfe=2^n2-(2^ia-2^(ia-1)+1)+1;

    % % Making Gaussian window
    % W=gausswin(floor((2*nx)/((2^(ia-1)+1))),2.5);
    
    % % Hanning window
    % W = hamming(floor((2*nx)/((2^(ia-1)+1))));  % window
    % W(1) = W(1)/2;  % repair constant-overlap-add for R=(M-1)/2
    % W(floor((2*nx)/((2^(ia-1)+1)))) = W(floor((2*nx)/((2^(ia-1)+1))))/2;
    
    % % Making B-Spline window
    W=b_spline_window(floor((nx)/((2^(ia-1)))),2);
    
    % Padding zero for the rest of the window
    WW=zeros(1,nx);
    WW(1:floor((nx)/((2^(ia-1)))))=W;
    FWW=fftshift(fft(WW));
    % Choosing only the middle part of the Fourier window
    WFWW=FWW(nx/2-2^(ia-2)+1:nx/2+2^(ia-2));
    WFWW=WFWW/max(abs(WFWW));
    
    % Just using a box car as a frequency window
    %WFWW=ones(1,length(WFWW));
    
    % Positive frequencies
    K(pfs:pfe)=(1/sqrt(2^(ia-1)))*fft((K(pfs:pfe).*abs(WFWW)));

    % Negative frequencies
    K(nfs:nfe)=(1/sqrt(2^(ia-1)))*fft((K(nfs:nfe).*abs(WFWW)));

end
