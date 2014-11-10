function yi=pwlint(x,y,xi)
% PWLINT: piecewise linear interpolation (much faster than interp1)
%
% yi=pwlint(x,y,xi)
%
% PWLINT performs linear interpolation when it is known that that
% input function (x,y) is piecewise linear. If length(x) is much less than
% the length(xi), this is MUCH faster than the built in INTERP1. Points
% in xi which are outside the bounds of x will return nans.
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

if(length(x)<=length(xi))
    nsegs=length(x)-1;
    yi=nan*zeros(size(xi));

    for k=1:nsegs

        %find the points in this line segment
        ii=between(x(k),x(k+1),xi,2);

        if( ii )
            % interpolate
            yi(ii)=y(k)+(y(k+1)-y(k))*(xi(ii)-x(k))/(x(k+1)-x(k));
        end

    end
else
    yi=nan*zeros(size(xi));
    for k=1:length(xi)
        ii=surround(x,xi(k));
        if(~isempty(ii))
            yi(k) = y(ii)*(x(ii+1)-xi(k))/(x(ii+1)-x(ii))+y(ii+1)*(xi(k)-x(ii))/(x(ii+1)-x(ii));
        else
            yi(k)=nan;
        end
    end
end
