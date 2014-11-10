function namenewloginit(transfer,q1,v1,q2,v2,q3,v3,q4,namelist,...
			zstart,zend,topsflag)
% namenewloginit(transfer,q1,v1,q2,q3,v3,q4,namelist,zstart,zend,topsflag)
%
% Initiate a dialog to name a new log and associate it with a well.
%   transfer= a string matlab command to be evaluated with EVAL at the 
%             completion of the dialog. Usually this will re-invoke the 
%             calling program
%	q1 = string giving the question to be asked to provide the name
%	v1 = the default value of the name. May be ''
%	q2 = string with the question asking for the name of the new well
%   	v2 = default value of the name of the new well. May be ''
%	q3 = question asking the inline coordinate of the new well
%	v3 = default value of the inline coordinate. May be ''
%   	q4 = fourth question associated with the namelist
% namelist = vector of names. If present, a toggle button will be present
%            allowing the new log to be associated with an existing well. If 
%            the namefield and inline coordinate fields will disappear to 
%            be replaced by a popup menu allowing the selection of one of 
%            the names in the vector.
%            Namelist should be in the format returned by 
%            objget(object,'fieldnames')
%   zstart = start depth of the log
%     zend = end depth of the log
% topsflag = if 1, then the tops are to be imported
%            if 0, then not imported
%            if -1, then there are no tops to worry about
%
% G.F. Margrave November 1993
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
	[m,n]=size(namelist);
	if( m > 1 )
		error(' namelist must be a vector ');
	end
	
	l=max([length(q1) length(q2) length(v1) length(transfer) length(v2) ...
			length(q3) length(v3) length(q4) n]);
	
	qs=ones(10,l);
	
	qs(1,1:length(q1))=q1;
	
	if(~isempty(v1))
		qs(2,1:length(v1))=v1;
	end
	
	qs(3,1:length(q2))=q2;
	
	qs(4,1:length(transfer))=transfer;
	
	if(~isempty(v2))
		qs(5,1:length(v2))=setstr(v2);
	end
	
	qs(6,(1:length(q3))) = q3;
	if(~isempty(v3))
		qs(7,(1:length(v3))) = v3;
	end
	qs(8,(1:length(q4))) = q4;
	if(~isempty(namelist))
	qs(9,1:n)=namelist;
	end
	
	qs(10,1)=zstart;
	qs(10,2)=zend;
	qs(10,3)=topsflag;
	
	set(gca,'userdata',qs);
	
	namenewlog('init');
        
