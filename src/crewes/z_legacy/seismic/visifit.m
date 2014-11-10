function [p,coeff]= visifit(y,x,norder)
% [p,coeff]= visifit(y,x,norder)
%
% VISIFIT allows the interactive fitting of an n'th order 
% polynomial to the data y(x). After plotting y(x) VISIFIT 
% graphical input mode and waits for the user to click both
% ends of the x coordinate range to be used in the fit.
%
% G.F. Margrave
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
x=x(:)';
y=y(:)';
% plot y(x)
 plot(x,y,'.')
 grid;
 text(.1,.05,'Select x coordinate range','sc');
% get the domain of the fit
[x1,y1]=ginput;
 while ~isempty(x1)
  xmin=min(x1);
  xmax=max(x1);
  if xmin==xmax, error(' Zero length domain chosen'),end
% fit the polynomial
  infit=near(x,xmin,xmax);
  coeff=polyfit(x(infit),y(infit),norder);
  p=polyval(coeff,x);
% draw new graph
   ind = find( p>1.1*max(y) );
   p(ind) = ones(size(p(ind)))*max(y);
   ind2 = find( p<min(y) );
   p(ind2) = ones(size(p(ind2)))*min(y);
  %plot(x,[y;p])
   plot(x,y,'.',x,p,'m')
   hold
   plot(x(ind),p(ind),'*m');
   plot(x(ind2),p(ind2),'*m');
  grid;
  text(.1,.00,'Select x coordinate range','sc');
  [x1,y1]=ginput;
end
