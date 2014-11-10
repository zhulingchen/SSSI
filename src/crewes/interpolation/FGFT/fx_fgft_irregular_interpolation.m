function Y=fx_fgft_irregular_interpolation(X,H,nh)
% This function interpolates regularly sampled seismic traces using
% a fast generalized Fourier transform. The Spitz method of using the low
% frequency portion for the interpolation of high frequencies is used.
%
%% Input
%   X: Data in t-x domain (with zero in the location of missing traces.
%
%
%% Output
%   Y: Interpolated data
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

[nt,nx]=size(X);

% to find the nearest power of two to padd with zeros
NT=2^(nextpow2(nt));
NX=2^(nextpow2(nh));

FX1=fft(X,NT,1);

% Starting loop over frequencies to be reconstructed
FY=zeros(NT,NX);
for ia=1:(NT/2)+1
    %Selecting frequency slices 
    x1=FX1(ia,:);
    
    ia 
    iter_cg=2;
    iter_bl=5;  
    
    % For irregular sampling
    INTD=irls_fitting_fgft(x1,H,NX,iter_cg,iter_bl);
    
    % Inverting D2 using inverse fgft and saving it in final array
    FY(ia,:)=INTD;
end

% Making the symmetric part of the frequency domain
FY(NT/2+2:NT,:)=conj(flipud(FY(2:NT/2,:)));

% Inverse Fourier transform from f-x to t-x domain
INTD=real(ifft(FY,NT,1));
Y=INTD(1:nt,1:nh);

% Normalizing trace values
for is=1:nh
    Y(:,is)= (Y(:,is) - mean(Y(:,is)))./std(Y(:,is));
end
