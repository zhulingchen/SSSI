function logout=logblock(login,zlog,zbdys,flag)
% logout=logblock(login,zlog,zbdys,flag)
%
% Block a well log.
%
% login = column vector containing the log samples
% zlog = column vector containing the sample depths
% zbdys = vector of block boundaries. There must be at least
%	2 boundaries. Log will be blocked between each pair of
%	boundaries. Ends of log will not be blocked unless boundaries
%   are defined at the ends. Boundaries will be sorted into ascending
%	order so they need not be on input
% flag = 1 ... MEAN : each log segment will be replaced by the log
%              mean value.
%        2 ... MEDIAN : each log segment will be replaced by the log
%              median value.
%        3 ... LINEAR TREND : each log segment will be replaced by 
%              a linear trend determined by least squares
% Note that in all cases NAN's are ignored unless the entire log
% segment is NAN in which case NAN is returned.
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
if(length(zbdys)<2)
	error('Must have at least 2 defined boundaries');
end
zbdys=sort(zbdys);
logout=login;
nblks=length(zbdys)-1;
for k=1:nblks
	iblk=between(zbdys(k),zbdys(k+1),zlog,1);
	%iblk indexes the points in the zone. Now find the subset
	%of them that are live
	ilive=find(~isnan(login(iblk)));
	if(isempty(ilive))
		logout(iblk)=nan*ones(size(iblk));
	else
		lseg=login(iblk(ilive));
		zseg=zlog(iblk(ilive));
		%determine mean, median or trend
		if(flag==1)
			p=zeros(1,2);
			p(2)=mean(lseg);
		elseif(flag==2)
			p=zeros(1,2);
			p(2)=median(lseg);
		elseif(flag==3)
			if(length(zseg)>=2)
				p=polyfit(zseg,lseg,1); 
			else
				p=[lseg(1) 0];
			end
		end
		logout(iblk)=polyval(p,zlog(iblk));
	end
end
