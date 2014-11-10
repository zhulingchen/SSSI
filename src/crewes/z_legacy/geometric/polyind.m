function [jmin,jmax] = polyind(ipoly,x,y)
% [jmin,jmax] = polyind(ipoly,x,y)
%
% x,y	= input arrays of polygon coordinates
% ipoly = input desired index of polygon is 
% 
% POLYIND returns:
% jmin 	= index of the first point
% jmax 	= index of the last point
%   
% Tom Bishop, Oct.93
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
% compute distmax
distmax=.05;    %threshold of 5 percent of picture
scalex=1./(max(x)-min(x));
scaley=1./(max(y)-min(y));
j=0;
distnext=distmax;
%
%loop thru all polygons
%
for loop = 1:ipoly
  jmin = j+1;  %reset from last polygon
  dist = distmax + 1;   % dont want 2nd point 
  j=jmin+1;
  while (dist > distmax) | (dist >= distnext) 
    j=j+1;
    if(j > length(x))
%      fprintf('reached end of array, length = %d\n',length(x))
      jmax = length(x);
      return
    end
    distx = scalex*(x(j)-x(jmin));
    disty = scaley*(y(j)-y(jmin));
    dist  = sqrt(distx.^2 + disty.^2);
    jnext = min(j+1,length(x));  %dont run off end of array
    distxnext = scalex*(x(jnext)-x(jmin));
    distynext = scaley*(y(jnext)-y(jmin));
    distnext  = sqrt(distxnext.^2 + distynext.^2);
%    fprintf('dist,distnext,j =  %8.4f %8.4f %d\n',...
%             dist,distnext,j) 
  end
  jmax=j;
end
