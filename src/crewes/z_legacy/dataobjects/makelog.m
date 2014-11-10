function wellog=makelog(name,data,flag,metric,logtype)

% wellog=makelog(name,data,flag,metric)
%
% Make a random Earth Object containing a well log.
% The storage conventions used here are documented in logobj.txt
%
%	name =  string naming the log
%	data = m rows by 2 columns. First column is either the depth or
%		time coordinate of a log sample and the second column is
%		the actual log samples.
%	flag = 0 ... its in depth
%	       1 ... its in time
%	metric = 0 ... the log is in imperial units
%	metric = 1 the log is in metric units
%	logtype = integer flag giving the log type as documented in
%		LAS2LOGTYPE or logobj.txt
%		default is -1 (unknown log type)
%
% by G.F. Margrave
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

if(nargin<5)
	logtype=-1;
end

wellog=objrand(name,'samples',data(:,2));
if(flag)
	wellog=objset(wellog,'time',data(:,1));
	wellog=objset(wellog,'datatype','tlog');
else
	wellog=objset(wellog,'depth',data(:,1));
	wellog=objset(wellog,'datatype','zlog');
end

%set the metric flag
wellog=objset(wellog,'dely',metric);

%set logtype
wellog=objset(wellog,'xnot',logtype);
