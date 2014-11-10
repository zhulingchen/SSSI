function clearpicks(figno)
% CLEARPICKS: clear (delete) the picks in a figure
%
% clearpicks(figno)
%
% CLEARPICKS deletes all picks from the current figure and removes
% them from the global PICKS. See also CLEARRAYS and CLEARLINES
%
% figno ... the figure number to delete from
% ********** default = gcf ***********
%
%
% by G.F. Margrave, CREWES, June 2000
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
if(nargin<1) figno=gcf; end

global PICKS

%resolve figure number
[nfigs,nnn]=size(PICKS);
doit=0;
for kkk=1:nfigs
    thisfig=PICKS{kkk,1};
    if(thisfig==figno)
        doit=1;
        break
    end
end
if(~doit)
    error('invalid figure number')
end

handles=PICKS{kkk,3};

hax=get(figno,'currentaxes');
hkids=get(hax,'children');
for k=1:length(hkids)
	if(strcmp(get(hkids(k),'type'),'line'))
% 		x=get(hkids(k),'xdata');
% 		if(length(x)==2)
% 			np=size(PICKS,1);
% 			y=get(hkids(k),'ydata');
% 			p=ones(np,1)*[x(1) y(1)];
% 			test=sum((PICKS-p)');
% 			ind=find(test==0);
% 			if(~isempty(ind))
% 				delete(hkids(k));
% 				PICKS(ind,:)=[];
% 			else
% 				p=ones(np,1)*[x(2) y(2)];
% 				test=sum((PICKS-p)');
% 				ind=find(test==0);
% 				if(~isempty(ind))
% 					delete(hkids(k));
% 					PICKS(ind,:)=[];
% 				end
% 			end
% 		end
        ind=find(handles==hkids(k));
        if(~isempty(ind))
            delete(hkids(k))
        end
	end
end

PICKS{kkk,2}=[];
PICKS{kkk,3}=[];
