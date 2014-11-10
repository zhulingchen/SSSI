function INTD=ls_mask_fitting_fgft(D,H,M,nh,iter_cg)
% This function fits the known mask function and available data using fgft.
%
%  Author(s): Mostafa Naghizadeh (University of Calgary)
%  Copyright 2010 Naghizadeh
%
%%% Note: Variables with upper-case letters are matrices or vector and
%%%       variables with lower-case letters are scalars
%
%%% Input:
% D: 1D Data (just available samples)
% Vector containig available data (offsets)
% M: Mask function in generalized Fourier domain
% nh: total number of offsets
% iter_cg: Number of Conjugate Gradient steps
%
%
%%% Output:
% INTD: Interpolated data
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

X= zeros(1,nh);
S=[D];
R=forward_operator(S,H,M,nh);
P=R;
rr_new = R*R';
err_old=S*S';
% Internal loop for CG fitting
for j = 1:iter_cg;
    Q=adjoint_operator(P,H,M,nh);
    alpha=rr_new/(Q*Q');
    X=X+alpha*P;
    S=S-alpha*Q;
    R=forward_operator(S,H,M,nh);
    rr_old=rr_new;
    rr_new=R*R';
    beta=rr_new/rr_old;
    P=R+beta*P;
    % error term computations
    err_new=S*S';
    err_old = err_new ;
end
INTD= sqrt(nh)*inverse_fgft(X);

%%%%%%%%%#################%%%%%%%%%%%%
function R=forward_operator(S,H,M,nh)
% Adjoint operator for band_limited reconstruction
TMP1=zeros(1,nh);
TMP1(H)= S;
R = (1/(sqrt(nh)))*forward_fgft(TMP1).*M; %length of nh


%%%%%%%%%#################%%%%%%%%%%%%
function Q=adjoint_operator(P,H,M,nh)
% Forward operator for band-limited reconstruction
TMP1 = (sqrt(nh))*inverse_fgft(M.*P);
Q = TMP1(H); %length ns
