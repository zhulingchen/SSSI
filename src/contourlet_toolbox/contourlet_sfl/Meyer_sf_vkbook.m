%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Meyer_sf_vkbook.m
%	
%   First Created: 08-26-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = Meyer_sf_vkbook(x)

%   The smooth passband function for constructing Meyer filters
%   See the following book for details:
%   M. Vetterli and J. Kovacevic, Wavelets and Subband Coding, Prentice
%   Hall, 1995.
  
theta = 3 * x .^ 2 - 2 * x .^ 3;
theta(x <= 0) = 0;
theta(x >= 1) = 1;