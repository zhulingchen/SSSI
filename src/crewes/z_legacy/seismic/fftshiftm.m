function y = fftshiftm(x)
% y= fftshiftm(x)
%
% This corrects the version of FFTSHIFT supplied in the Matlab toolbox
% which is erroneous for matricies. For an input matrix [1 2;3 4],
% FFTSHIFT returns [4 3;1 2] while FFTSHIFTM returns [3 4;1 2]
%
% Correction by G.F. Margrave, May 1991
%
%FFTSHIFTM Shift FFT.  For vectors FFTSHIFT(X) returns a vector with the
%	left and right halves swapped.  For matrices, FFTSHIFT(X) swaps
%	the first and third quadrants and the second and fourth quadrants.
%	FFTSHIFT is useful for FFT processing, moving the zeroth lag to
%	the center of the spectrum.
%
%	J.N. Little 6-23-86
%	Copyright (c) 1986 by the MathWorks, Inc.
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
[m,n] = size(x);
m1 = 1:floor(m/2);
m2 = floor(m/2)+1:m;
% Note: n/2+1 references the Nyquist. so the spectrum after FFTSHIFT 
% goes from -Nyquist -> DC -> one sample before + Nyquist
n1 = 1:floor(n/2);
n2 = floor(n/2+1):n;
% Note: can remove the first two cases when null handling is fixed.
if m == 1
	y = [x(n2) x(n1)];
elseif n == 1
	y = [x(m2); x(m1)];
else
 y = [x(m2,n1) x(m2,n2); x(m1,n1) x(m1,n2)];
end
