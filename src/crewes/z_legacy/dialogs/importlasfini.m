function [lognums,lognames,wellname,x,zstart,zend,topsflag,units,kb]...
          =importlasfini
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
%test for cancel
a=get(gca,'userdata');
if(a==-1)
    lognums=-1;
    lognames=-1;wellname=-1;x=-1;zstart=-1;
    zend=-1;topsflag=-1;units=-1;kb=-1;
    return;
end
ind=find(isnan(a));
lognums=abs(a(1:ind(1)-1));
lognames=setstr(a(ind(1)+3:ind(2)-1));
m=abs(a(ind(1)+1));
n=abs(a(ind(1)+2));
lognames=reshape(lognames,m,n);
wellname=setstr( char( a(ind(2)+1:ind(3)-1)) );
units=setstr(a(ind(3)+1:ind(4)-1));
x= abs(a(ind(4)+1));
zstart=abs(a(ind(4)+2));
zend=abs(a(ind(4)+3));
topsflag=abs(a(ind(4)+4));
kb=abs(a(ind(4)+5));
