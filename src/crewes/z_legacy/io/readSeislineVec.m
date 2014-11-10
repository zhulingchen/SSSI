function obj = readSeislineVec(filename)
%   obj = readSeislineVec(filename)
%   builds Gary's object from seisline vector card file
%   card file is named filename, it is built with extract_vector
%   see also: fold
%
%  T. N. BISHOP,  OCTOBER 1993,  CPTC CANADA
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
fid = fopen(filename,'r');
%A = fscanf(fid,'%*s %d %f %f %f %f');
A = fscanf(fid,'%*s %f %f %f %f %f');
B = reshape(A,5,length(A)/5)';
trno = B(:,1);
trno = -trno;
trno = trno-min(trno);
x = B(:,2);
y = B(:,3);
s = sqrt((x-x(1)).^2 + (y-y(1)).^2);
t =B(:,4);
t = -t;
lmax = 0;
fprintf ('\n min and max dist = %10.1f %10.1f \n',min(s),max(s));
fprintf ('min and max time = %10.1f %10.1f \n',min(t),max(t));
vel = input('enter velocity, use same units as dist, time, above');
t = t*vel;     %convert to depth
fprintf ('min and max depths = %10.1f %10.1f \n',min(t),max(t));
mx = 0;
npoly = 0;
while mx < length(t)  % read in polys till read end of arrays
  npoly = npoly + 1;
  [mn,mx]=polyind(npoly,s,t);
  l = mx-mn+1;
  if(l > lmax)  lmax = l; end 
end
lmax = 1.1 * lmax;
fprintf ('number of polygons = %5d \n',npoly);
for ipoly=1:npoly
  [mn,mx]=polyind(ipoly,s,t);
  l = mx-mn+1;
  column1 = [s(mn:mx);NaN*ones(lmax-l,1)];
  column2 = [t(mn:mx);NaN*ones(lmax-l,1)];
  if(ipoly == 1) 
    obj=gridobj('Fold',column1,[int2str(ipoly),'s']);
    obj=objset(obj,[int2str(ipoly),'t'],column2);
  else
    obj=objset(obj,[int2str(ipoly),'s'],column1);
    obj=objset(obj,[int2str(ipoly),'t'],column2);
  end
end
