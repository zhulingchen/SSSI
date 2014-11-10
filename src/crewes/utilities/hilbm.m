function y = hilbm(x)
% HILBM: Hilbert transform
%
% HILBERT Hilbert transform.
%
% Modified from Matlab's version to require the input to be
% a power of 2 in length. by GFM 
%
%	HILBERT(X) is the Hilbert transform of the real part
%	of vector X.  The real part of the result is the original
%	real data; the imaginary part is the actual Hilbert
%	transform.  See also FFT and IFFT.

%	Charles R. Denham, January 7, 1988.
%	Revised by LS, 11-19-88.
%	Copyright (C) 1988 the MathWorks, Inc.

% Reference: Jon Claerbout, Introduction to
%            Geophysical Data Analysis.
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

test=2.^nextpow2(length(x));
 if test~= length(x)
  error(' input vector length must be a power of 2')
 end
yy = fft(real(x));
m = length(yy);
if m ~= 1
	h = [1; 2*ones(m/2,1); zeros(m-m/2-1,1)];
	yy(:) = yy(:).*h;
end
y = ifft(yy);
