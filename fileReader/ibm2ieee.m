function d = ibm2ieee (ibmf)

% Name: ibm2ieee
% Abstract: convert a matrix of IBM/360 32-bit floats
% to IEEE doubles.
%
% IBMF is the matrix of IBM/360 32-bit floats each
% stored as a 32 bit unsigned big-endian integer
% in a MATLAB double.
%
% The format of a IBM/360 32-bit float is:
%
% sign 7-bit exponent 24 bit fraction
% The base is 16. The decimal point is to
% the left of the fraction. The exponent is
% biased by +64.
%
% The basic idea is that you use floating point on
% the various fields.
%
% ieee = sign * 16 ^ (exponent - 64) * fraction / 16 ^ 6
%
% By: Martin Knapp-Cordes
% MathWorks, Inc.
%
% Date(s): Jun 95 - 28, 29

% $Revision: 1.2 $ $Date: 1995/06/29 14:50:03 $
%----------------------------------------------------------------------------
%
if (nargin ~= 1)
error ('Wrong number of arguments.');
elseif (isempty(ibmf))
error ('Argument is an empty matrix');
end
%
aibmf = sprintf('%08x',ibmf);
%
% hexd(1) - 1st hex digit - first bit is sign, next 3 bits high order exponent
% hexd(2) - 2nd hex digit - bits of low order exponent
% hexd(3) - 3rd-8th hex digit - bits of fraction
%
hexd = sscanf(aibmf,'%1x%1x%6x',[3,inf]);
d = (1 - (hexd(1,:) >= 8) .* 2) .* ...
16 .^ ((hexd(1,:) - (hexd(1,:) >= 8) .* 8) .* 16 + hexd(2,:) ...
- 70).* hexd(3,:);
d = reshape(d,size(ibmf));