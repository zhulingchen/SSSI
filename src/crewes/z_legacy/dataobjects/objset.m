function object=objset(objin,item,value,arg4)

% object = objset( objin, item, value);
% object = objset( objin, item, value, datanum);
% object = objset( objin, item, value, dataname);
%
% OBJSET accepts an EarthObject (created by RANDOBJ, GRIDOBJ)
% and sets the values of one of it's items. If the
% item already exists, then its values are replaced; otherwise
% a new item is added to the object and its values inserted.
% New items are assumed to be data (as opposed to headers) and must be
% the same size and geometry as other data items already in the object. The
% identity of header fields in the object is fixed so new header information
% cannot be added; however, the existing headers can be changed (updated) at
% any point. A data item may be protected from accidental overwrite by setting its
% protect flag to 'on'. The first form of the command is used most often
% and sets items that are not specific to a particular dataset within the 
% object. Examples: (let myObj be an existing object)
%    myObj = objset(myObj, 'username', 'barney'); ... change the username to barney
%    yourObj = objset(myObj,'username','barney'); ... as above but create a new object
%    myObj = objset(myObj,' Leduc ',newGrid); ... a grid referred to by the matlab
%		variable "newGrid" is put in the object and named ' Leduc '
%    myObj = objset(myObj,'protect','on',' Leduc '); ... the ' Leduc ' grid is 
%		protected.
%    myObj = objset(myObj,'protect','off',' Leduc '); ... the ' Leduc ' grid is 
%		un-protected.
%    myObj = objset(myObj,'protect','on',3); ... the third grid is un-protected.
%
% objin = the variable name of the object to be updated
% item = a string specifying the item to be set (max of 30 characters)
% value = the value(s) of the item
% datanum = the sequential number of one of the data fields in the object
% dataname = string giving the name of one of the data fields in the object
%
% object = the output object. Normally this will be the same matrix as objin
%
% For more information type: help earthobj
%
%
% by G.F. Margrave, November 1993
% Updated July 1997
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

% determine which type of object this is

type = oget(objin,'objtype');

% branch accordingly

if( strcmp(type,'cont') )
	if(nargin == 4)
		object=ocset(objin,item,value,arg4);
	else
		object=ocset(objin,item,value);
	end
else
	if(nargin == 4)
		object=oset(objin,item,value,arg4);
	else
		object=oset(objin,item,value);
	end
end
