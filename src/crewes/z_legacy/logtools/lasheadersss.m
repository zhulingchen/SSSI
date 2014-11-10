function lasheader = lasheadersss(lash,start,stop,step,units)

% lasheader = lasheadersss(lasheader,start,stop,step,units)
% lasheader = lasheadersss(lasheader,start,stop,step)
%
% Given an las header in a string matrix (such as is provided by
% readlas) LASHEADERSSS change its start, stop, and step to
% the appropriate values. If units is provided, it is used to specify
% the units of start stop and step, otherwise they are unchanged.
%
% G.F. Margrave, CCR, Aug 94
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
    units='';
end
[m,n]=size(lash);

for k=1:m
    tmpstr=lash(k,:);
    [t,r]=strtok(tmpstr);
    if(findstr(t,'STRT'))
        lash(k,:)=replacenum(tmpstr,start,units);
    end
    if(findstr(t,'STOP'))
        lash(k,:)=replacenum(tmpstr,stop,units);
    end
    if(findstr(t,'STEP'))
        lash(k,:)=replacenum(tmpstr,step,units);
    end
end
lasheader=lash;

function tmpstr=replacenum(instring,num,units)
        n=length(instring);
        [t,r]=strtok(instring);
        %test for units adjacent to the flag
        [tt,rr]=strtok(t,'.');
        if(isempty(rr))
            [tt,r]=strtok(r); %find units
            iu=findstr(tt,instring);
            if(~isempty(units))%replace units
                instring=[instring(1:iu(1)-1) '.' units instring(iu(1)+2:end)];
            end
            numfield=strtok(r);
        else
            iu=findstr(rr,instring);
            if(~isempty(units))%replace units
                instring=[instring(1:iu(1)-1) '.' units instring(iu(1)+2:end)];
            end
            numfield=strtok(r);
        end
        if(isempty(numfield)|strcmp(numfield(1),':'))
            %if here, then "numfield" is actually a comment and we have no start specified
            % ugh
            tmpstr=[instring(1:iu(1)+2) num2str(num) instring(iu(1)+3:end)];
        else
            ii=findstr(numfield,instring);%find the numeric field
            tmpstr=[instring(1:ii(1)-1) num2str(num) instring(ii(1)+length(numfield):end)];
            
        end
        if(length(tmpstr)<n) tmpstr=[tmpstr blanks(n-length(tmpstr))]; end
        if(length(tmpstr)>n) tmpstr=tmpstr(1:n); end
