function W=b_spline_window(L,m)
%
% This function creates a B-Spline window by consecutive convolution of box car
% functions. The order 1 B-Spline window is a triangular one. Order 0 is a
% box car function. This function is used to create windows so that in the case
% of half overlap between their lobes the partition of unity is satisfied.
%
%%% Input:
%   L= Length of window (Better to be power of 2)
%   m= Order of the B-Spline window (Integer value)
%
%%% Output
%   W=B-Spline window normalized to 1.
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

% Starting window size. Also order 0 B-Spline window
W=ones(1,L);

BW=W;
if (m>0)
    for ia=1:m
        % Convolving the latest order of B-Spline window with box car function
        TW=conv(W,BW);
        % Normalizing the window
        TW=TW/max(TW);
        % Getting new B-Spline window by decimating the convolution result
        W=TW(1:2:2*L);
    end
end
