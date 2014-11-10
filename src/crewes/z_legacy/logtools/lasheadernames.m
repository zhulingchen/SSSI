function lasheader = lasheadernames(lasheader,names,flag,desc)
%
% lasheader = lasheadernames(lasheader,badnames)
% lasheader = lasheadernames(lasheader,goodnames,1)
%
% Given an las header in a string matrix (such as is provided by
% readlas) LASHEADERNAME searches it for references to curves whose
% 4 letter mnemonics are in badnames and remove all such references.
% Or if a third argument whose value is 1 is supplied, then the names
% are assumed to be added to the header.
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
%
%
if(nargin<3)
	flag=0;
end
if(~flag)
  badnames=names;
  k=1;
  str=lasheader(k,:);
  [nrows,ncols]=size(lasheader);
  [nbad,m]=size(badnames);
	  %find the curve block
  while( ~strcmp(lower(str(1:4)),'~cur') )
		k=k+1;
		str=lasheader(k,:);
	 end
%
	 %step through the curve block deleting any lines with badnames
	 k=k+1;
	 str=lasheader(k,:);
	 goodnames=[];
	 while(str(1,1)~='~')
		kbad=0;
		if(str(1)==' ') ichk=2:5; else ichk=1:4; end
		for kk=1:nbad
			if(strcmp(str(ichk),badnames(kk,:)))
				kbad=kk;
				break;
			end
		end
		if(kbad)
			lasheader(k,:)=[];
		else
			k=k+1;
			if(str(1)~='#')
				goodnames=[goodnames '  ' str(ichk)];
			end
		end
		str=lasheader(k,:);
	 end
	  [nlines,nchar]=size(lasheader);
	  %tidy up the last entry
	  str=upper(lasheader(nlines,:));
	  cols=length(str);
	  %ind=findstr(str,'DEPT');
	  %if(isempty(ind))
	%	ind=findstr(str,'TIME');
	%  end
	%  if(isempty(ind))
	%	ind=findstr(str,'ETIM');
	%  end
	%  if(isempty(ind))
	%	error('logic failure in lasheadernames');
	%  end
	%  i2= find(str(ind(1):cols)==' ');
	  str=['~A  ' goodnames];
		lasheader(nlines,:)=[str blanks(cols-length(str))];
	else
		goodnames=names;
		k=1;
		str=lasheader(k,:);
		[nrows,ncols]=size(lasheader);
		[ngood,m]=size(goodnames);
		if(nargin<4)
			desc=setstr(32*ones(ngood,1));
		end
		[nn,ndesc]=size(desc);
		%find the curve block
		while( ~strcmp(lower(str(1:4)),'~cur') )
			k=k+1;
			str=lasheader(k,:);
		end
		k=k+1;
		str=lasheader(k,:);
		%step to the end of the block
		%count curves
		ncurves=0;
		if(str(1,1)==' ')
			ichk=2:5;
		else
			ichk=1:4;
		end
		metric=-1;%assume time log
		while( ~strcmp(str(1),'~') )
			%determine depth units
			if(strcmp(str(ichk),'DEPT'))
				if(str(ichk(4)+2)=='M')
					metric=1;%metric log
				else
					metric=0;%imperial log
				end
			end
			if(~strcmp(str(1),'#'))
				ncurves=ncurves+1;
			end
			k=k+1;
			str=lasheader(k,:);
		end
		%copy the last curve record as a template
		k=k-1;
		str=lasheader(k,:);
		ind=find(str==':');
		if(~isempty(ind))
			str(11:ind(1)-1)=blanks(ind(1)-11);
		else
			ind=27;
			str(11:ind(1)-1)=blanks(ind(1)-11);
		end
		str(ind+1:length(str))=blanks(length(str)-ind);
    if(str(1)==' ') ifill=2:5; else ifill=1:4; end
		imnem=ifill(4)+2:ifill(4)+5;
		%loop over good names adding one line per each
		for kk=1:ngood
			str(ifill)=goodnames(kk,:);
			str(imnem)=lasmnem2units(goodnames(kk,:),metric);
			ncurves=ncurves+1;
			ns=int2str(ncurves);
			bb=blanks(3-length(ns));
			str(ind+6:ind+5+length(desc(kk,:))+3)=[ns bb desc(kk,:)];
			lasheader=[lasheader(1:k,:); str; lasheader(k+1:nrows,:)];
			k=k+1;
			nrows=nrows+1;
		end
		
		%add the names to the last string
		str=lasheader(nrows,:);
		ind=find(abs(str)~=32);
		imax=max(ind);
		for kk=1:ngood
			str(imax+1:imax+5)=[' ' goodnames(kk,:)];
			imax=imax+5;
		end
		
		lasheader(nrows,:)=str(1:ncols);
			
	end
