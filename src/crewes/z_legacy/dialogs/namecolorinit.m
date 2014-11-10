function namecolorinit(transfer,q1,v1,q2,v2,q3,namelist)
% Please use "namecolorle" when developing new programs.
%
% namecolorinit(transfer,q1,v1,q2,v2,q3,namelist)
% namecolorinit(transfer,q1,v1,q2,v2)
%
% Initiate a name/color selection dialog.
% transfer = a string matlab command to be evaluated with EVAL at the 
%            completion of the dialog. Usually this will re-invoke the 
%            calling program
%	q1 = string giving the question to be asked to provide the name
%	v1 = the default value of the name. May be ''
%	q2 = string giving the question asking for a color.
%       v2 = default value of the color as an rgb triplet. [.5 .5 .5] 
%            gives a neutral grey
%       q3 = optional third question associated with the namelist
% namelist = optional vector of names. If present, a toggle button will 
%            be present allowing the current object to be named the same 
%            as a previous one. If clicked, the namefield and color 
%            selection will disappear to be replaced by a popup menu
%            allowing the selection of one of the names in the vector. 
%            Namelist should be in the format returned by 
%            objget(object,'fieldnames') 
%
% G.F. Margrave November 1993
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
if( nargin == 6 )
   error('incorrect number of arguments');
end
if(nargin<6)
   q3=[];
   namelist=[];
end
if(nargin<5)
   error('minimum of 5 arguments required');
end
[m,n]=size(namelist);
if( m > 1 )
   error(' namelist must be a vector ');
end
l=max([length(q1) length(q2) length(v1) length(transfer) length(v2) ...
	length(q3) n]);
if( m == 0)
   qs=ones(5,l);
else
   qs=ones(7,l);
end
qs(1,1:length(q1))=q1;
qs(2,1:length(v1))=v1;
qs(3,1:length(q2))=q2;
qs(4,1:length(transfer))=transfer;
qs(5,1:length(v2))=v2;
if( m>0 )
   qs(6,(1:length(q3))) = q3;
   qs(7,1:n)=namelist;
end
set(gca,'userdata',qs);
namecolor('init');
