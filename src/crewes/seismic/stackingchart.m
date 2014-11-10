function stackingchart(xshot,xrec,dx)
% STACKINGCHART ... plots seismic line geometry
% 
% stackingchart(xshot,xrec)
%
% xshot ... vector of shot coordinates
% xrec ... cell array of receiver coordinates
%   Requirement: length(xshot) must equal length(xrec)
% dx ... cmp binning interval
%  **** default half of the median receiver spacing *****
% 
%
%
% G.F. Margrave, CREWES Project, August 2013
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
if(length(xshot)~=length(xrec))
    error('xshot and xrec must be the same length')
end

if(~iscell(xrec))
    error('xrec must be a cell array');
end

if(iscell(xshot))
    error('xshot must not be a cell array')
end
xcmpmin=10^20;
xcmpmax=-xcmpmin;
if(nargin<3)
    dxs=zeros(size(xshot));
    for k=1:length(xshot)
        dxs(k)=mean(diff(xrec{k}));
        xcmps=.5*(xshot(k)+xrec{k});
        if(min(xcmps)<xcmpmin)
            xcmpmin=min(xcmps);
        end
        if(max(xcmps)>xcmpmax)
            xcmpmax=max(xcmps);
        end
    end
    dx=mean(dxs)/2;
end
xcmp=xcmpmin:dx:xcmpmax;
figure
nshots=length(xshot);
for k=1:nshots
    h=line(xrec{k},xshot(k)*ones(size(xrec{k})));
    set(h,'linestyle','none','marker','.','color','r');
    h=line(xshot(k),xshot(k));
    set(h,'linestyle','none','marker','*');
end
xlabel('Receiver coordinate')
ylabel('Shot coordinate')
title('Shot-Receiver Chart')

figure
fold=zeros(size(xcmp));
for k=1:nshots
    xmid=(xshot(k)+xrec{k})/2;
    xoff=abs(xrec{k}-xshot(k));
    h=line(xmid,xoff);
    set(h,'linestyle','none','marker','x');
%     h=line(xshot(k),xshot(k));
%     set(h,'linestyle','none','marker','*');
    %count fold
    icmp=ceil((xmid-xcmpmin)/dx);
    for kk=1:length(icmp)
        if(icmp(kk)>0 && icmp(kk)<=length(fold))
            fold(icmp(kk))=fold(icmp(kk))+1;
        end
    end
end
xlabel('Midpoint coordinate')
ylabel('Offset coordinate')
title('Midpoint-Offset Chart')

figure
plot(xcmp,fold)
title('CMP Fold')
xlabel('cmp coordinate')