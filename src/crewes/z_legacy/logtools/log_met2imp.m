function logout=log_met2imp(login,logtype)
% logout=log_met2imp(login,logtype)
%
% LOG_MET2IMP converts a wellog from metric to imperial units.
% The intent is to recognize all log types and convert their units
% correctly. If you feel this is not being properly done, please
% contact the programer/author listed below.
%
%	login ... input log object
%	May be an earth object as described in logobj.txt or a 2 column matrix
%	where the first column is depth and the second the log samples
% 	logtype ... one of the following strings: 'sonic' 'density'
%       or a numeric code as described in logobj.txt. The best way is to use
%	the numeric code because that allows the specification of lots of
%	different log types. The two ttext strings are maintained as a 
%	nod to the past.
%
%	logout ... output log object or 2 column matrix. Will be the same
%		kind of object as was input
%
% G.F. Margrave, June 1994
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
logout=login;
if(~isstr(logtype))
	if(logtype==0) logtype='sonic';
	elseif( logtype==1 | logtype==2 | logtype==3) logtype='density';
	elseif( logtype==6 )
		logtype='caliper';
	elseif( logtype==7 )%s and p sonics are the same for unit changes
		logtype='sonic';
	elseif( logtype==13 | logtype==14 |logtype==15 )
		logtype='resistivity';
	else
		logtype='whatever';
	end
end
objin=0;
if( isearthobj(login) )
	objin=1;
	%get the metric flag
	flag=objget(login,'dely');
	if(flag==0)
		error(' Log is already imperial');
	end
	logout=objset(logout,'dely',1);
	% get the samples
	samps=objget(login,'samples');
	z=objget(login,'depth');
	t=objget(login,'time');
else
	z=login(:,1);
	samps=login(:,2);
end
if(~isempty(z))
		z=z/(.3048);
end
if( strcmp(lower(logtype),'sonic') )
	samps=samps/3.280840;
elseif( strcmp(lower(logtype),'density'))
	samps=samps/1000;
elseif( strcmp(logtype,'caliper') )
	samps=samps/25.4;
elseif( strcmp(logtype,'resistivity') )
	samps=samps/.3048;
end
if( objin )
		logout=objset(logout,'samples',samps);
		if(~isempty(z))
			logout=objset(logout,'depth',z);
		end
	else
		logout=[z(:) samps(:)];
end
