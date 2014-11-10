function Y = matrix_scaling(X, nx2, ny2)
% This function projects the values of an [nx1 by ny1] array into an array
% of size [nx2 by ny2] using the nearest neighbor interpolation.
%
%% Input
%   X: The original array
%   nx1: number of samples in x direction for Output matrix
%   ny1: number of samples in y direction for Output matrix
%
%% Output
%   Y: Output matrix
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


% Determining the size of the original data
[nx1,ny1]=size(X);

% Vectors to determine the x,y coordinates of the new matrix related to old grid
XC=round((nx1/nx2)*[1:nx2]);
YC=round((ny1/ny2)*[1:ny2]);

% Make sure indices of matrix are matched
for ic=1:nx2
    if (XC(ic)==0)
        XC(ic)=1;
    end
    if (XC(ic)>nx1)
        XC(ic)=nx1;
    end
end
for ic=1:ny2
    if (YC(ic)==0)
        YC(ic)=1;
    end
    if (YC(ic)>ny1)
        YC(ic)=ny1;
    end
end

% Replacing values of X into Y
Y=X(XC,YC);
