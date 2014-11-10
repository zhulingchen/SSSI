function ind = closetraceazi(xxx,yyy,xw,yw,rd,azi) 
% function ind = closetraceazi(xxx,yyy,xw,yw,rd,azi) 
%
%  finds closest trace (point) in line to given well location (xw,yw)
%  (as long as the trace is within a radius rd)
%  the distance MUST be defined along the azimuth azi, given in degrees
%  line defined by points xxx(i),yyy(i)
%  returns index of closest trace 
%  returns null if no trace is close enough
%
%  T. N. BISHOP,  DECEMBER 1993,  CPTC CANADA
%   see also seis2well, linewelltie
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
    ind = [];
    dclose= closeline(xxx,yyy,xw,yw);
    if dclose <= rd
      mm = tan(2*pi*azi/360.);
      bb = yw - mm*xw;
      dmin = inf;
      m = diff(yyy)./diff(xxx);
      b = yyy(2:length(yyy)) - m.*xxx(2:length(xxx));
      xc = (bb - b)./(m - mm);
      yc = m.*xc + b;
      for i = 1:(length(xxx)-1) 
        i1 = i + 1;
        if ~isempty(between(xxx(i),xxx(i1),xc(i),2)) 
          if ~isempty(between(yyy(i),yyy(i1),yc(i),2)) 
%   find closest trace (point in xxx,yyy) to point
            if abs(xc(i)-xxx(i)) < abs(xc(i)-xxx(i1))
 	      ind = i;
            else
	      ind = i1;
 	    end
            d=sqrt( (xxx(ind)-xw).^2 + (yyy(ind)-yw).^2 );
            if(d < dmin)
              dmin = d;
      	      indmin = ind;
	    end
          end
        end
      end
%      plot(xxx(indmin),yyy(indmin),'c*')  plot closest trace point
      if dmin <= rd   %min.dist of all points is less than radius
        ind = indmin;
      else      %min.dist along azimuth not less than radius
        ind = [];
      end
    else
      ind = [];   %not close enough to even check
    end
