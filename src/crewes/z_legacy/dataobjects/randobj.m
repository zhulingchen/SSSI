function E = randobj(objname,data,dataname,x,y,datatype)

% E= randobj( name, data, dataname, x, y, datatype)
% E= randobj( name, data, dataname, x, y)
% E= randobj( name, data, dataname, x)
% E= randobj( name, data, dataname)
% E= randobj( name, data)
%
% RANDOBJ creates an EarthObject designed to store 
% earth science data sampled at arbitrary (x,y) locations. 
% The returned matrix, E, should NEVER be used directly in 
% computations because it is designed as a general data storage 
% object. Instead, use the appropriate OBJGET and OBJSET calls to 
% get and set the items in E.
%
% name = string vector specifying the object's name 
% data = vector containing the data
% x = vector containing the x coordinates of the data
%	 ******* default= no field set *********
%
% y = vector containing the y coordinates of the data
%        ******* default= no field set *********
%
% datatype = string vector specifying the data type
% 	 ******* default = '    ' *********
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
     error('You must supply at least two arguments');
 end
% set defaults
 if nargin < 6 
     datatype = '    ';
 end
 if nargin<5
     y = [];
 end
if nargin < 4
        x=[];
end
 if nargin<3
    datatype=data;
    data=[];
    dataname=[];
 end  

% make sure data is a vector
if( min(size(data)) > 1 )
	error(' Data must be a vector not a matrix');
end

% declare E
 E= buildobj(length(data),length(objname));
 
[m,n]=size(data);

% set object version number
 E.objvernum=1.0;
%set object creation date
 E.objcreate=fix(clock);
% set user field
 E.username='matlab';
%set x axis name
 E.xaxisname='x axis';
%set y axis name
 E.yaxisname='y axis';
%set azimuth of north
 E.north=90.0;
% set Azimuth of east
 E.east=0.0;
% set global origin
 E.globorigin=[0.,0.];
% set scale factor
 E.xyscale=1.0;
% set rotation angle
 E.xyrotate=0.;
% set flag for randobj
 E.objtype='rand';
% set data type
 E.datatype=datatype;
% set datarows
 E.datarows=m;
% set datacols
 E.datacols=n;
% set object's name
 E.name=objname;
 
% set the object's x coords
if( length(x)>0 )
  if( length(x)==max(m,n) )
    E=setfield(E,'x',x);
  else
    error(' incorrect data size for x coordinates');
  end
end

% set the y coords
if( length(y)>0 )
  if( length(y)==max(m,n) )
    E=setfield(E,'y',y);
  else
    error(' incorrect data size for y coordinates');
  end
end
 
% set the data name
 E.data=setfield(E.data,'dataname',dataname);
% set data create 
 E.data=setfield(E.data,'datacreate',fix(clock));
% set data modified
 E.data=setfield(E.data,'datamodified',fix(clock));
% set the data
if( length(data)>0 )
 E.data=setfield(E.data,dataname,data);
end

