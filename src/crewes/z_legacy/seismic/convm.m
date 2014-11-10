function c = convm(a, b)
% c= convm(a,b)
%
% convm is a modification of the 'conv' routine from the MATLAB
% toolbox. The changes make it more convenient for seismic purposes
% in that the output vector, c, has a length equal to the first
% input vector,  a. Thus, 'a' might correspond to a reflectivity
% estimate to be convolved with a wavelet contained in 'b' to
% produce a synthetic seismic response 'c'. It is assumed that
% the wavelet in b is causal and that the first sample occurs at time zero.
% For non-causal wavelets, use 'convz'. An abort will occur if
% b is longer than a. CONVM will correctly handle an ensemble matrix.
% 
%CONV	Convolution and polynomial multiplication.
%	C = CONV(A, B) convolves vectors A and B.  The resulting
%	vector is length LENGTH(A)
%	If A and B are vectors of polynomial coefficients, convolving
%	them is equivalent to multiplying the two polynomials.
%	See also XCORR and DECONV and CONVZ
%
%	J.N. Little 4-21-85
%	Revised 9-3-87 JNL
%	Copyright (c) 1985, 1987 by the MathWorks, Inc.
%
% modified to 'convm' by G.F. Margrave, May 1991
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
[rows,cols]=size(a);  
if rows>cols, a=a.';end
na = max(size(a));
nvecs= min(size(a));
nb = max(size(b));
% Convolution, polynomial multiplication, and FIR digital
% filtering are all the same operations.  Since FILTER
% is a fast built-in primitive, we'll use it for CONV.
% This routine assumes that the first argument is the 
% which is longer than the wavelet (second argument). An
% abort occurs if this is not so.
%
if na<nb, error(' First vector must be longer than second'),end
if nvecs==1,
 
	  if na > 1
		  a(na) = 0;
	  end 
  	c = filter(b, 1, a);
   mw=mwhalf(na,4);
   c=c.*mw;
 else
  c=zeros(size(a));
  mw=mwhalf(na,4); 
  for k=1:nvecs
    c(k,:)=filter(b,1,a(k,:)).*mw;
  end
 end
if rows>cols, c=c.';end
