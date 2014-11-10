function [lmin,lmax,n] = split(line)
%
% [lmin,lmax,n] = split(line)
%
% line = input vector of character variable, usually the line name
%
% lmin	= vector of the indicies of the min.trace in each line
% lmax = vector of the indicies of the max.trace in each line
% n 	= number of lines
%  
% for example, if line = [aa aa aa bb bb cc cc cc], then
% lmin = [1 4 6], lmax = [3 5 8], n=3
%
% T. N. BISHOP,  NOVEMBER 1993,  CPTC CANADA
%   
% SEE ALSO SEIS2WELL, LINEWELLTIE, PLOTSEISWELL
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
  i = 1;
  n = 1;
  lmin(n) = 1;
  lmax(n) = 1;
  imax=length(line);
while (n < 100 )
  while( line(i,:) == line(i+1,:) ) %while no change in line name
    i=i+1;
    lmax(n) = i;
    if (i+1) > imax  %end of file
      return
    end
  end
  i = i+1;
  n = n+1;
  lmin(n)=lmax(n-1)+1;
end
disp('you have more than 100 lines');
