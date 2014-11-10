function lasheader_t = lasheadertime(lasheader,logt,ids)
%
% lasheader_t = lasheadertime(lasheader,logt,ids)
%
% Given an las header in a string matrix (such as is provided by
% readlas) LASHEADERTIME converts its start, stop, and interval to
% the appropriate times. Also, it is assumed that the logs to be 
% represented by the header are a subset of those originally 
% described by the header in depth. The argument 'ids' gives those
% the logs numbers (e.g. 1,2 etc) which are to be represented and
% this program deletes any other log id information
%
% T.N. Bishop, CCR, July 94
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
  lasheader_t=lasheader;
  [nlines,nchar]=size(lasheader);
%
%  need to modify the lasheader, STRT, STOP, STEP
%       also DEPT? DEPTH?
%
  k=1;
	 nlogs=length(ids);
	 flag=0;
  str=lasheader_t(k,:);
  [nrows,ncols]=size(lasheader_t);
  badnames=zeros(30,4);
  nbad=0;
  n=0;
  while( ~strcmp(str(1:2),'~A') )
%
	  i4=2:5;
	  if(str(1)~=' ')
		i4=1:4;
	  end
	  if(strcmp(str(i4),'STRT'))
	    ind=find(str==':')-1;
	    startt= min(logt(:,1));   %min time value
	    ttt=num2str(startt);
	  end
	  if(strcmp(str(i4),'STOP'))
	    ind=find(str==':')-1;
	    stopt= max(logt(:,1));   %max time value
	    ttt=num2str(stopt);
	  end
	  if(strcmp(str(i4),'STEP'))
	    ind=find(str==':')-1;
	    stept= logt(2,1)-logt(1,1);   %step time value
	    ttt=num2str(stept);
	  end
	  if(strcmp(str(i4),'STRT')| ...
		strcmp(str(i4),'STOP') | ...
		strcmp(str(i4),'STEP') )
	    bl='                                          ';
	    ttt = ['  ' ttt bl(1:ind-11-length(ttt)) ];
	    str(10:ind)=ttt;
	    str(max(i4)+2)='S';
	  end
	  if(flag)
		  if(str(1) ~= '~')
			n=n+1;
			ind=find(ids==n);
			if(isempty(ind))
				nbad=nbad+1;
				badnames(nbad,:)=str(i4);
				str=[];
			else
				it=find(str==':');
				if(~isempty(it))
					it2=find(abs(str(it:length(str)))==32);
					tmp=[num2str(ind+1) str(it+it2(2)-1:length(str))];
					str=[str(1:it+it2(1)-1) tmp];
				end
			end
		  else
			flag=0;
		  end
	  end
	  if(isempty(str))
		  if(strcmp(str(i4),'DEPT'))
			 str(i4(1):max(i4)+3)='ETIM.S ';
			 ind=findstr(str,'DEPTH');
			 if(~isempty(ind))
				str(ind:ind+4)='TIME ';
			 end
			 % also the next n lines id the logs
			 flag=1;
			 n=0;
		  end
		end
		
		if(length(str)>=ncols)
   	  lasheader_t(k,:) = setstr(str(1:ncols));
		elseif( length(str)>0 )
			lasheader_t(k,:) = setstr([str 32*ones(1,ncols-length(str))]);
		else
			lasheader_t(k,:)=str;
		end
		  if(~isempty(str))
			  k=k+1;
			end
		  str=lasheader_t(k,:);
  end
  %tidy up the last entry
  cols=length(str);
  for kbad=1:nbad
		ind=findstr(str,badnames(kbad,:));
		str=[str(1:ind-1) str(ind+11:length(str))];
	 end
	 ind=findstr(str,'DEPTH');
	 if(~isempty(ind))
		str(ind:ind+4)='TIME ';
	 end
	 if( length(str)<cols)
		 lasheader_t(k,:)=setstr([str 32*ones(1,cols-length(str))]);
	 else
		 lasheader_t(k,:)=str(1:cols);
	 end
