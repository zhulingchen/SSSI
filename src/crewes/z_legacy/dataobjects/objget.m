function itemValue = objget( object, item, arg3)

% itemValue = objget(object, item)
% itemvalue = objget(object, item, datanum)
% itemvalue = objget(object, item, dataname)
%
% OBJGET performs the inverse operation to OBJSET in that it
% retrieves data from an EarthObject (created by RANDOBJ, GRIDOBJ)
% Most invocations will use the first syntax. Examples:
%	objget(myobject,'name') ... returns the object's name
%	objget(myobject,'x') ... returns the object's x coordinates
%   objget(myobject,'fred') ... returns the data field named fred
%	objget(myobject,5) ... returns the fifth dataset in the object
%
% The second syntax is provided to obtain fields specific to of a dataset 
% when only its name or number are known. For example, the name of the 
% 5th dataset is obtained with:
%   objget(myobject,'dataname',5) ... will return the name of the 5th data grid
% If the returned name is 'peace river' for example, then its last 
% modification date can be obtained by either of:
%   objget(myobject,'datamodified',5); ... 'peace river' is the 5th dataset
%   objget(myobject,'datamodified','peace river'); ... ask for it my name
%
% object = the EarthObject to be interrogated
% item = a string identifying the item to be retrieved, or an integer
%          denoting the index of the item. This last case is useful if 
%          a query such as objget( object, 'fieldnames') has first been
%          made to determine the names of available fields in the object.
%          Any field can then be obtained by simply specifying it's number
%          (the first is number 1 etc...)
% datanum = an integer referring to one of the data fields in the object. 
%           Valid numbers are 1 through m where m is the number of data 
%           fields stored.  
%           (and the number of rows returned by objget(myObject,'namesmatrix'))
% dataname = a string giving the name of a dataset stored in the object. 
%            Blanks are important: 
%            the string 'leduc' will not match ' leduc' or 'leduc '
%
% itemValue =  the returned value of the item. The form of itemValue
%     may be any of scalar, vector, or matrix depending on what
%     was stored.
%
% by G.F. Margrave, November 1993.
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

% determine the type of object and branch

type=object.objtype;

if( strcmp(type,'cont') )
	if( nargin==3 )
		itemValue = ocget(object, item, arg3);
	else
		itemValue = ocget(object, item);
	end
else
	if( nargin==3 )
		itemValue = oget(object, item, arg3);
	else
		itemValue = oget(object, item);
	end
end
