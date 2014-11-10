function D=scale_fgft(D)
% This function removes the sqrt(nh) factor from forward transformed FGFT
% data
%
%% Input
%   D: 1D array that contains the coefficients of FGFT (1 x nk)
%
%% Output
%   D: Scales data
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
nk=length(D);

%PK=zeros(nk,nk);

n2=nextpow2(nk);

for ia=2:n2-1
    % At each ia step we will take the fft of positive and
    % negative frequencies with variant window size.
    pfs=2^ia-2^(ia-1)+2;
    pfe=2^ia+1;
    nfs=(2^n2-(2^ia-2^(ia-1)+1)+1)-2^(ia-1)+2;
    nfe=2^n2-(2^ia-2^(ia-1)+1)+2;
    
    
    % Positive frequencies
    D(pfs:pfe)=sqrt(2^(ia-1))*D(pfs:pfe);
        
    % Negative frequencies
    D(nfs:nfe)=sqrt(2^(ia-1))*D(nfs:nfe);
    
end
