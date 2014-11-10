function E = objrand(objname,dataname,data,x,y,datatype)

% E= objrand( name, dataname, data, x, y, datatype)
% E= objrand( name, dataname, data, x, y)
% E= objrand( name, dataname, data, x)
% E= objrand( name, dataname, data)
% E= objrand( name, datatype)
%
%
% OBJRAND creates an EarthObject designed to store 
% earth science data sampled at arbitrary (x,y) locations. 
% The returned matrix, E, should NEVER be used directly in 
% computations because it is designed as a general data storage 
% object. Instead, use the appropriate OBJGET and OBJSET calls to 
% get and set the items in E.
%
% name = string vector specifying the object's name 
%	    requirement: length(name) <= length(data)
% data = vector containing the data
% x = vector containing the x coordinates of the data
%				 ******* default= no field set *********
% y = vector containing the y coordinates of the data
%				 ******* default= no field set *********
% datatype = string vector specifying the data type
%                    requirement: length(datatype) == 4 
%     ******* default = '    ' *********
%
% E =  the returned EarthObject
%
% by G.F. Margrave, March 1993
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
     error('You must supply at least two arguments');
 end
% set defaults
 if nargin < 6 
     datatype = '    ';
 end
 if nargin<5
     y = [];
 end
 if nargin<4
     x = [];
 end  
if nargin < 3
 datatype=dataname;
 data=[];
 dataname=[];
end

% make sure data is a vector
if( min(size(data)) > 1 )
	error(' Data must be a vector not a matrix');
end
if( ~ischar(dataname) & ~ischar(datatype) )
	error(' second argument must be a string ');
end

E=randobj(objname,data,dataname,x,y,datatype);
