function a=askthingsfini
% Please use "askthingsle" if developing new programs
%
% a=askthingsfini
%
% ASKTHINGSFINI terminates the askthings dialog and returns a matrix of
% strings containing the answers to the asked questions. One per row. If the
% user canceled then a is equal to -1.
%
% G.F. Margrave January 1994
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
as=get(gca,'userdata');
% hfig=abs(as(length(as))); %fix
% as(length(as))=[]; %fix
% test for cancel
if(as==-1)
	a=as;
	return;
end
ind=find(double(as)==255);
% determine the maximum answer length
test=[0 ind];
len=diff(test)-1;
maxlen=max(len);
% get the answers
na=length(ind);
a=setstr(32*ones(na,maxlen));
kbegin=1;
for k=1:na
	a(k,1:len(k))=as(kbegin:ind(k)-1);
	kbegin=ind(k)+1;
end
set(gca,'userdata',[]);
% close(hfig); %fix
