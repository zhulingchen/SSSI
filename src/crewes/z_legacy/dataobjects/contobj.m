function E = contobj(objname,subobj,subobjname,datatype,subobjdatatype)

% E= contobj( name, subobj,subobjname,datatype,subobjdatatype)
% E= contobj( name, subobj,subobjname,subobjdatatype)
% E= contobj( name, subobj,subobjname)
% E= contobj( name, datatype)
%
%
% CONTOBJ creates a container EarthObject designed to store other 
% EarthObjects. The last syntactical form creates an empty container object.
% The returned matrix, E, should NEVER be used directly in 
% computations because it is designed as a general data storage 
% object. Instead, use the appropriate OBJGET and OBJSET calls to 
% get and set the items in E.
%
% name = string vector specifying the object's name 
%	    requirement: length(name) <= length(data)
% subobj = another object to be contained in this one
% subobjname = the name of the subobj
%     ************ default = objget(subobj,'name') ************
%
% datatype = string vector specifying the data type of the container object
%     ******* default = '    ' *********
%
% subobjdatatype = string vector specifying the data type of the subobject.
%     ********* length(subobjdatatype) = 4 ******************** 
%     ********* default = objget(subobject,'datatype') unless the
%	  subobject is a simple scalar, matrix or vector. In that case the
%         default depends on size(subobject)
%		size(subobject)		default
%		    [1,1]		 'sclr'
%		    [n,1]		 'cvec'
%		    [1,n]		 'rvec'
%		    [m,n]		 'mtrx'
%
% E =  the returned EarthObject
%
% by G.F. Margrave, March 1993
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

if nargin < 2
     error('You must supply at least three arguments');
end

% declare E
 E= buildcontobj;
 
% set defaults
if nargin < 4 
     datatype = '    ';
end
if nargin < 3
   datatype = subobj;
   subobj = [];
   subobjname = [];
end 
[m,n]=size(subobj);
if nargin < 5
	if( isearthobj(subobj) )
	          E.data.dataobjtype = objget(subobj,'objtype');
 	elseif( m==1 & n==1)
 		  E.data.dataobjtype = 'sclr';
 	elseif( n>1 & m==1 )
 		  E.data.dataobjtype = 'rvec';
 	elseif( n==1 & m>1 )
 		  E.data.dataobjtype = 'cvec';
 	elseif( n>1 & m>1 )
 		  E.data.dataobjtype = 'mtrx';
 	end
 end
 
% set object version number
 E.objvernum=1.0;
%set object creation date
 E.objcreate=fix(clock);
% fix object modification date
E.objmodified=fix(clock);
% set user field
 E.username='matlab';
% set flag for contobj
 E.objtype='cont';
% set data type
 if( length(datatype)~=4 )
 	error(' data types must be 4 characters in length');
 else
 	E.datatype=datatype;
 end
% set object's name
 E.name=objname;

% set data modified values
 E.data.datamodified=fix(clock);
  
% set data name
 E.data.dataname=subobjname;
 
% set protection value (default=0)
 E.data.protect=0;
 
% set datadatatype
if( ischar(subobj) )
 datatype='strg';
else
 datatype='smpl';
end
 E.data.datadatatype=datatype;
   
% set the data
if( ~isempty(subobj) )
 	E.data=setfield(E.data,subobjname,subobj);
end
