function PK=plot_fgft(K)
% This function plots the coefficients of a generalized Fourier transform. 
% It transforms 1D coefficients into a 2D array. 
%
%% Input
%   K: 1D array that contains the coefficients of FGFT (1 x nk)
%
%% Output
%   PK: 2D array (nk x nk)
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


% Determining the size of coefficients
nk=length(K);

%PK=zeros(nk,nk);

n2=nextpow2(nk);

% Upscaling first and second coefficients
PKp1=matrix_scaling(K(1),1,nk);
PKp2=matrix_scaling(K(2),1,nk);
PKn1=matrix_scaling(K(nk-1),1,nk);
PKn2=matrix_scaling(K(nk),1,nk);

% Starting loop to pick parts of K and upscale them
PPK=[];
NPK=[];
for ia=2:n2-1
    % At each ia step we will take the fft of positive and
    % negative frequencies with variant window size.
    pfs=2^ia-2^(ia-1)+1;
    pfe=2^ia;
    nfs=(2^n2-(2^ia-2^(ia-1)+1)+1)-2^(ia-1)+1;
    nfe=2^n2-(2^ia-2^(ia-1)+1)+1;
    
    
    % Positive frequencies
    PPK=[PPK; matrix_scaling(K(pfs:pfe),2^(ia-1),nk)];
        
    % Negative frequencies
    NPK=[matrix_scaling(K(nfs:nfe),2^(ia-1),nk); NPK];
    
end

PK=fliplr([PKp1;PKp2;PPK;NPK;PKn1;PKn2]);
