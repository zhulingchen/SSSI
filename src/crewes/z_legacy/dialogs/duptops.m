function [ret1,ret2]=duptops(arg1,arg2,arg3)
% duptops(tops,topnames,transfer) will initiate the dialog
% [tops,topnames]=duptops('fini') returns the culled tops
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
 
if(~isstr(arg1))
	action='init';
	tops=arg1;
	topnames=arg2;
	transfer=arg3;
	arg1=[];
	arg2=[];
else
	action=arg1;
end
if(strcmp(action,'init'))
 %search for duplicate tops
 [tops,ind]=sort(tops);
 ntops=length(tops);
 topnames=topnames(ind,:);
 dtops=diff(tops);
 ind=find(dtops==0.0);
 ndup=0;
 ret1=1;
 if( isempty(ind) )
  	 ind=0;
	 ret1=-1;
	 return;
 end
 %let idup be the index of the starting element of each set of dups
 %let ndup be the number of duplicates in each set
 %idup=[];
 %ndup=[];
 idup=ind;
 ndup=zeros(size(ind));
 ikill=ndup;
 ireal=1;
 ikeep=ones(size(tops));
	
	for k=1:length(ind)
		test=tops-tops(ind(k));
		itest=find(test==0);
		ndup(k)=length(itest);
		%test to see if we've already counted this one
		if(k>1)
			if(idup(k)<=idup(ireal)+ndup(ireal)-1)
				ikill(k)=1;
			else
				%ireal=ireal+1;
				ireal=k;
				ik=zeros(1,ndup(k));
				ik(1)=1;
				ikeep(idup(k):idup(k)+ndup(k)-1)=ik;
			end
		else
			ik=zeros(1,ndup(k));
			ik(1)=1;
			ikeep(idup(k):idup(k)+ndup(k)-1)=ik;
		end
	end
	if( ikill > 0 )
		idup(ikill)=[];
		ndup(ikill)=[];
	end
	
	hfig=figure('visible','off');
	pos=get(hfig,'position');
	figwid=400;
	fight=100;
	set(hfig,'position',[pos(1:2) figwid fight],'visible','on');
	%done and cancel buttons
	sep=.01;
	xnow=.01;
	ynow=.01;
	width=.15;
	height=.2;
	hdone=uicontrol('style','pushbutton','string','Done',...
		'units','normalized','position',[xnow ynow width height],...
		'callback','duptops(''done'')','foregroundcolor','r',...
		'userdata',ikeep);
	xnow=xnow+width+sep;
	hcancel=uicontrol('style','pushbutton','string','Cancel',...
		'units','normalized','position',[xnow ynow width height],...
		'callback','duptops(''cancel'')','userdata',transfer);
	msg1='On the left are depths which have duplicate tops.';
	msg2='For each depth select the ONE top to keep on the right.';
 xnow=sep;
 ynow=ynow+5*sep+height;
 height=1*height;
 width=1-2*sep;
	hmsg2=uicontrol('style','text','string',msg2,'max',1,...
		'units','normalized','position',[xnow ynow width height],...
		'userdata',ndup);
 ynow=ynow+height;
	hmsg1=uicontrol('style','text','string',msg1,'max',1,...
		'units','normalized','position',[xnow ynow width height],...
		'userdata',idup);
	ynow=ynow+height+5*sep;
	xnow=2*sep;
	width=.5-3*sep;
	%lbl=[];
	lbl=setstr(32*ones(length(idup),20));
	for k=1:length(idup)
		%lbl=strmat(lbl,num2str(tops(idup(k))));
		tmp=sprintf('%6.4f',tops(idup(k)));
		lbl(k,1:length(tmp))=tmp;
	end
	htops=uicontrol('style','popupmenu','string',lbl,...
		'units','normalized','position',[xnow ynow width height],...
		'callback','duptops(''chgz'')','backgroundcolor','c',...
		'userdata',tops);
	xnow=xnow+width+2*sep;
	htopnames=uicontrol('style','popupmenu',...
		'string',topnames(idup(1):idup(1)+ndup(1)-1,:),...
		'units','normalized','position',[xnow ynow width height],...
		'callback','duptops(''choose'')','backgroundcolor','c',...
		'userdata',topnames);
	set(hfig,'visible','on');
	set(hfig,'userdata',[hdone,hcancel,hmsg1,hmsg2,htops,htopnames]);
	%
	% userdata assignments
	% hdone ... ikeep
	% hcancel ... transfer function
	% hmsg1 ... idup
	% hmsg2 ... ndup
	% htops ... tops
	% htopnames ... topnames
 return;
end
if(strcmp(action,'chgz'))
	h=get(gcf,'userdata');
	hdone=h(1);
	hmsg1=h(3);
	hmsg2=h(4);
	htops=h(5);
	htopnames=h(6);
	ikeep=get(hdone,'userdata');
	idup=get(hmsg1,'userdata');
	ndup=get(hmsg2,'userdata');
	tops=get(htops,'userdata');
	topnames=get(htopnames,'userdata');
	%determine the new depth
	inew=get(htops,'value');
	znew=tops(idup(inew));
	%get the names for this depth
	names=topnames(idup(inew):idup(inew)+ndup(inew)-1,:);
	%determine which one is currently kept
	ikeepnow=find(ikeep(idup(inew):idup(inew)+ndup(inew)-1));
	%set the names
	set(htopnames,'string',names,'value',ikeepnow);
	return;
end
if(strcmp(action,'choose'))
	h=get(gcf,'userdata');
	hdone=h(1);
	hmsg1=h(3);
	hmsg2=h(4);
	htops=h(5);
	htopnames=h(6);
	ikeep=get(hdone,'userdata');
	idup=get(hmsg1,'userdata');
	ndup=get(hmsg2,'userdata');
	tops=get(htops,'userdata');
	topnames=get(htopnames,'userdata');
	%get the value 
	ikeepnew=get(htopnames,'value');
	%save it in ikeep
	itop=get(htops,'value');
	temp=zeros(1,ndup(itop));
	temp(ikeepnew)=1;
	ikeep(idup(itop):idup(itop)+ndup(itop)-1)=temp;
	set(hdone,'userdata',ikeep);
	return;
end
%process the done button
if(strcmp(action,'done'))
	h=get(gcf,'userdata');
	hcancel=h(2);
	transfer=get(hcancel,'userdata');
	%call the transfer function
	eval(transfer);
	return;
end
%process the cancel button
if(strcmp(action,'cancel'))
	h=get(gcf,'userdata');
	hdone=h(1);
	hcancel=h(2);
	transfer=get(hcancel,'userdata');
	%set the cancel flag
	set(hdone,'userdata',-1);
	%call the transfer function
	eval(transfer);
	return;
end
if(strcmp(action,'fini'))
	h=get(gcf,'userdata');
	hdone=h(1);
	hmsg1=h(3);
	hmsg2=h(4);
	htops=h(5);
	htopnames=h(6);
	ikeep=get(hdone,'userdata');
	idup=get(hmsg1,'userdata');
	ndup=get(hmsg2,'userdata');
	tops=get(htops,'userdata');
	topnames=get(htopnames,'userdata');
	%test for cancel
	if(ikeep==-1)
		ret1=-1;
		ret2=[];
	else
		%cull the tops and topnames
		ret1=tops(ikeep);
		ret2=topnames(ikeep,:);
	end
	close(gcf);
	return;
end
