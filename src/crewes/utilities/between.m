function indicies = between(x1,x2,testpts,flag)
% BETWEEN: logical test, finds sampls in vector between given bounds
%
% indicies = between(x1,x2,testpts,flag)
% indicies = between(x1,x2,testpts)
%
% returns the indicies of those points in the array testpts which lie between
% the points x1 and x2. If no testpts are found between x1 and x2 then a
% single scalar 0 (false) is returned. Flag determines the exact nature of
% the inclusion of the endpoints:
%    if flag == 0, then not endpoints are included
%    if flag == 1, then x1 is included. i.e. if a test point is precisely
%                  equal to x1,	it will be considered "between"
%    if flag == 2, then both x1 and x2 are included
%  ******** flag defaults to 0 ****
%
% Function works regardless of whether x1 < x2 or x2 < x1.
%
% by G.F. Margrave, May 1991
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

if nargin < 4, flag =0; end

if(length(x1)~=1 || length(x2)~=1)
    error('x1 and x2 must be scalars')
end

if flag == 0
    if( x1 < x2)
        indicies = find( (testpts > x1)&(testpts < x2) );
        if(isempty(indicies))indicies=0;end
        return;
    else
        indicies = find( (testpts > x2)&(testpts < x1) );
        if(isempty(indicies))indicies=0;end
        return;
    end
end
if flag == 1
    if( x1 < x2)
        indicies = find( (testpts >= x1)&(testpts < x2) );
        if(isempty(indicies))indicies=0;end
        return;
    else
        indicies = find( (testpts > x2)&(testpts <= x1) );
        if(isempty(indicies))indicies=0;end
        return;
    end
end
if flag == 2
    if( x1 < x2)
        indicies = find( (testpts >= x1)&(testpts <= x2) );
        if(isempty(indicies))indicies=0;end
        return;
    else
        indicies = find( (testpts >= x2)&(testpts <= x1) );
        if(isempty(indicies))indicies=0;end
        return;
    end
end
