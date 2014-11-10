function y = filttri(x,n)
% y = filttri(x,n)
%
% THIS IS BASICALLY A WEIGHTED MOVING AVERAGE OPERATOR
% x is the input vector; the no.of pts. of the vector
%    need not be specified
% n is the number of points of the triangle smoother, 
% y is the smoothed x, convolved with the
%    triangle of n pts.  Because the added edges at
%    the beginning and end are truncated, y will have
%    the same number of points as x.
%  T.N.Bishop, CCR, 4/94
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
%  compute the triangle function of n points
n1 = fix(n/2)+1;
box = ones(1,n1)/n1;
f1 = conv(box,box);
		%pad input data with end values
npad = fix(n/2);
onepad = ones(1,npad);
xpad = [x(1)*onepad x x(length(x))*onepad];
%  compute the center of the boxcar, assume n is odd
nc = fix(n/2) + 1;
%  use Gary's fct from seis.toolbox, to get same no.of
%  output points as input
y = convz(xpad,f1,nc);
		% unpad output data with end values
y = y((npad+1):(length(y)-npad));
