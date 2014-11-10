function [shot_nmo,xcmp,offsets,fold]=nmor_cmp(shot,t,x,xshot,vrmsmod,xv,tv,dxcmp,x0cmp,x1cmp,flag,smax)
% NMOR_CMP: Remove NMO and map to CMP for a shot record
%
% [shot_nmo,xcmp,offsets,fold]=nmor_cmp(shot,t,x,xshot,vrmsmod,xv,tv,dxcmp,x0cmp,x1cmp,flag,smax)
%
% shot ... input shot record
% t ... time coordinate for shot
% x ... receiver coordinate for shot
% xshot ... source coordinate for shot
% vrmsmod ... rms velocity model for shot
% xv ... x coordinate for vrmsmod
% tv ... time coordinate for vrmsmod
% dxcmp ... cmp increment
% x0cmp ... initial cmp coordinate
% x1cmp ... final cmp coordinate (not used for flag=0)
% flag ... 0 means output cmp range just spans the input range
%         1 means output cmp range extends from x0cmp to x1cmp and ouput
%         gather is padded with zero traces
% smax ... NMO stretch mute (percent). Muting occurs when stretch is
% greater than this value
%  ************ default 30 ****************
%
% shot_nmo ... nmo corrected shot with traces mapped to cmp axis
% xcmp ... cmp axis for shot_nmo
% offsets ... mean offset of each output trace
% fold ... fold for shot_nmo (NOTE: shot_nmo is not fold normalized)
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
if(nargin<12)
    smax=30;
end

if(length(x)~=size(shot,2))
    error('invalid x coordinate for shot gather')
end
if(length(t)~=size(shot,1))
    error('invalid t coordinate for shot gather')
end
if(length(xv)~=size(vrmsmod,2))
    error('invalid x coordinate for velocity model')
end
if(length(tv)~=size(vrmsmod,1))
    error('invalid t coordinate for velocity model')
end
tv=tv(:);%force column vector
xv=xv(:)';%force row vector

%raw cmps and offsets
xcmpraw=(x+xshot)/2;
xoff=abs(x-xshot);

%expand velocity model if need be to span input cmps
if(min(xv)>min(xcmpraw))||(max(xv)<max(xcmpraw))
    xv=[min(xcmpraw)-dxcmp xv max(xcmpraw)+dxcmp];
    vrmsmod=[vrmsmod(:,1) vrmsmod vrmsmod(:,end)];%repeat first and last traces
end
%expand velocity model in time
if(tv(end)<t(end))
    tv=[tv;t(end)];
    vrmsmod=[vrmsmod;vrmsmod(:,end)];
end
%test to see if we need velocity interpolation in time
test=sum(abs(diff(tv)));
if(test>100*eps)
    interp_t=1;
else
    interp_t=0;
end
%make cmp axis

if(flag==0)
    x1cmp=(round((max(xcmpraw)-x0cmp)/dx)+1)*dx;
end
icmp=round((xcmpraw-x0cmp)/dxcmp)+1;%cmp bin numbers for each trace
xcmpnom=x0cmp:dxcmp:x1cmp;%nominal cmp bin centers from origin to max of this shot
if(flag==0)
    xcmp=xcmpnom(min(icmp):max(icmp));
    i0=min(icmp);%first output bin number
else
    xcmp=xcmpnom;
    i0=1;
end
shot_nmo=zeros(length(t),length(xcmp));
fold=zeros(size(xcmp));
offsets=zeros(size(xcmp));

%loop over input traces
for k=1:length(x)
    %determine velocity for this trace
    ind=surround(xv,xcmpraw(k));
    vrms=(xcmpraw(k)-xv(ind(1)+1))*vrmsmod(:,ind(1))/(xv(ind(1))-xv(ind(1)+1))...
        +(xcmpraw(k)-xv(ind(1)))*vrmsmod(:,ind(1)+1)/(xv(ind(1)+1)-xv(ind(1)));
    if(interp_t)
       vrms=interp1(tv,vrms,t);
    end
    %remove moveout on current trace
    params=nan*ones(1,3);
    params(1)=smax;
    shot_nmo(:,icmp(k)-i0+1)=nmor(shot(:,k),t,xoff(k),vrms,1,params);
    %the following line calculates the mean offset of the traces falling in a bin
    offsets(icmp(k)-i0+1)=(offsets(icmp(k)-i0+1)*fold(icmp(k)-i0+1)+xoff(k))/(fold(icmp(k)-i0+1)+1);
    fold(icmp(k)-i0+1)=fold(icmp(k)-i0+1)+1;%increment fold
    
end
%normalize by fold
% for k=1:length(xcmp)
%     if(fold(k)>0)
%         shot_nmo(:,k)=shot_nmo(:,k)/fold(k);
%     end
% end
    
