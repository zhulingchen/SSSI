function [ap,ae,tp,xp]=picker(seis,t,x,te,xe,delt,flag)
% PICKER ... pick amplitudes on seismic events
%
% [ap,ae,tp,xp]=picker(seis,t,x,te,xe,delt,flag)
%
% PICKER extracts amplitudes from an event in a seismic matrix. The event
% is defined by a single-valued traveltime trajectory. The picks are made
% at interpolated times as necessary. The defined trajectory is interpreted
% as spatially continuous and times are interpolated at each trace location
% if not specified. So a linear event need only have its endpoints
% specified. Trajectories are not extrapolated, so if a pick is desired at
% every trace then the ends of the trajectory should exceed the bounds of
% the seismic. The times specified for the trajectory may bounce around as
% needed to follow some perceived event.
%
% SEE ALSO ipick, which provides an interactiv interface to picker.
%
% seis ... seismic matrix
% t ... time coordinate vector for seis (length(t)=size(seis,1))
% x ... x coordinate vector for seis (length(x)=seis(seis,2))
% te ... vector of times defining the event
% xe ... vector of x's defining the event (length(te)=length(xe)). The
%   entries of xe must be monotonic to ensure the event is single valued.
% delt ... halfwidth of fairway around the defined trajectory
%    NOTE: set delt to zero to force the picks to be taken at exactly the 
%       trajectory times. This makes methods 3 and 4 identical and is not
%       available for methods 5, 6, and 7.
% flag ... 1: pick max(abs) amplitude in fairway
%          2: pick closest peak of Hilbert envelope in fairway
%          3: pick closest peak in fairway
%          4: pick closest trough in fairway
%          5: pick closest + to - zero crossing in fairway
%          6: pick closest - to + zero crossing in fairway
%          7: pick closest zero crossing of either polarity
%          8: pick first breaks by max(sta(env)/lta(env))
%          9: pick first breaks by max Hilbert envelope in fairway
% 
% ap ... vector of picked amplitudes (for options 5,6,7 this is the
%           amplitude of the Hilbert envelope at the zero crossing)
% ae ... vector of event amplitudes at the same times as the picks
% NOTE: ae and ap are identical for methods 3 and 4
% tp ... vector of pick times
% xp ... vector of pick x coordinates
%
% NOTE: The various methods behave slightly differently if the criteria are
% not met in the fairway (usually happens if delt is too small). 
% For methods 5,6,&7, if no suitable zero crossing is found in the fairway  
% then the pick is assigned nan (not a number). This will cause the graphic 
% display to be blank for that pick. For methods 3 and 4, if no peak is found 
% then the pick is assigned to the maximum value in the fairway. For method 4, 
% the pick is assigned to the minimum if no trough is found.
%
% G.F. Margrave, 2007-10, CREWES
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

global THRESH TIMES FBAVE

[nt,nx]=size(seis);
if(length(t)~=nt)
    error('time vector is the wrong size')
end
if(length(x)~=nx)
    error('coordinate vector is wrong size')
end
if(length(te)~=length(xe))
    error('te and xe must be the same length')
end
if(prod(diff(xe))<0)
    error('xe must be monotonic')
end

dt=t(2)-t(1);
if(delt==0)
    method='exact';
    delt=5*dt;%a nominal fairway
else
    method='interp';
end

if(delt<dt)
    warning('fairway specified is less than the sample size, rounding up')
end

te=te(:);%force column vector

%determine the part of the trajectory that lies in the matrix
tmin=min(t);
tmax=max(t);
xmin=min(x);
xmax=max(x);
ind=find((te<tmin) | (te>tmax));
if(~isempty(ind))
    te(ind)=[];
    xe(ind)=[];
end
ind=find((xe<xmin) | (xe>xmax));
if(~isempty(ind))
    te(ind)=[];
    xe(ind)=[];
end
if(isempty(te))
    ap=[];
    ae=[];
    tp=[];
    xp=[];
    return
end

%slice the matrix
indx=between(min(xe),max(xe),x,2);
xp=x(indx);
te_int=pwlint(xe,te,x(indx));%interpolate times at each trace
te_int=te_int(:);%column vector

samps=slicemat(seis(:,indx),(te_int-t(1))/dt+1,ceil(delt/dt));%grab the fairway slice
[nt2,nx2]=size(samps);
nh=(nt2-1)/2;
trel=(-nh:nh)*dt;%relative time of slice
%interpolate to 10 times the samples. So, pick times will be to nearest
%neighbor but at one-tenth the sample rate
dt2=dt/10;
trel2=(trel(1):dt2:trel(end))';%column vector
samps2=zeros(length(trel2),nx2);
for k=1:nx2
    samps2(:,k)=spline(trel,samps(:,k),trel2);
end
switch method
    case {'exact'}
        switch flag
            case {1}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    tabs=round(te_int(k)/dt)*dt+trel2;
                    ae(k)=spline(tabs,samps2(:,k),te_int(k));
                    ap(k)=abs(ae(k));
                    tp(k)=te_int(k);
                end
            case {2}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    tabs=round(te_int(k)/dt)*dt+trel2;
                    env=abs(hilbert(samps2(:,k)));
                    ae(k)=spline(tabs,samps2(:,k),te_int(k));
                    ap(k)=spline(tabs,env,te_int(k));
                    tp(k)=te_int(k);
                end
            case {3,4}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    tabs=round(te_int(k)/dt)*dt+trel2;
                    ae(k)=spline(tabs,samps2(:,k),te_int(k));
                    ap(k)=ae(k);
                    tp(k)=te_int(k);
                end
            case {5,6,7,8,9}
                msgbox('Exact picking not possible for cases 5,6,7,8,9')
            otherwise
                error('unknown picking method')
        end
    case {'interp'}
        switch flag
            case {1}
                [ap,ind]=max(abs(samps2));%find max abs value on each trace
                ae=zeros(nx2,1);
                for k=1:nx2
                    ae(k)=samps2(ind(k),k); %extract the true amp at the location of max_abs
                end
                tp=te_int+trel2(ind);%true time of the pick
            case {2}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    env=abs(hilbert(samps2(:,k)));
                    iex=findex(env,1);
                    if(~isempty(iex))
                        it=near(trel2(iex),0);%find closest to the event time
                        ap(k)=env(iex(it(1)));
                        ae(k)=samps2(iex(it(1)),k);
                        tp(k)=te_int(k)+trel2(iex(it(1)));
                    else
                        [ap(k),ind]=max(env);
                        tp(k)=te_int(k)+trel2(ind(1));
                        ae(k)=samps2(ind(1),k);
                    end
                end
            case {3}
                ae=zeros(nx2,1);
                ap=ae;
                tp=ae;
                for k=1:nx2
                    iex=findex(samps2(:,k),1);%find peaks
                    if(~isempty(iex))
                        it=near(trel2(iex),0);%find closest to the event time
                        ap(k)=samps2(iex(it(1)),k);
                        tp(k)=te_int(k)+trel2(iex(it(1)));
                        ae(k)=ap(k);
                    else
                        [ap(k),ind]=max(samps2(:,k));
                        tp(k)=te_int(k)+trel2(ind(1));
                        ae(k)=ap(k);
                    end
                end
            case {4}
                ae=zeros(nx2,1);
                ap=ae;
                tp=ae;
                for k=1:nx2
                    iex=findex(samps2(:,k),-1);%find troughs
                    if(~isempty(iex))
                        it=near(trel2(iex),0);%find closest to the event time
                        ap(k)=samps2(iex(it(1)),k);
                        tp(k)=te_int(k)+trel2(iex(it(1)));
                        ae(k)=ap(k);
                    else
                        [ap(k),ind]=min(samps2(:,k));
                        tp(k)=te_int(k)+trel2(ind(1));
                        ae(k)=ap(k);
                    end
                end
            case {5}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    ind_xc=find_zero_crossings(samps2(:,k),1);
                    env=abs(hilbert(samps2(:,k)));
                    if(~isempty(ind_xc))
                        it=near(trel2(ind_xc),0);%find closest to the event time
                        tp(k)=te_int(k)+trel2(ind_xc(it(1)));
                        ae(k)=samps2(ind_xc(it),k);
                        ap(k)=env(ind_xc(it));
                    else
                        tp(k)=nan;
                        ae(k)=nan;
                        ap(k)=nan;
                    end
                end
            case {6}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    ind_xc=find_zero_crossings(samps2(:,k),-1);
                    env=abs(hilbert(samps2(:,k)));
                    if(~isempty(ind_xc))
                        it=near(trel2(ind_xc),0);%find closest to the event time
                        tp(k)=te_int(k)+trel2(ind_xc(it(1)));
                        ae(k)=samps2(ind_xc(it),k);
                        ap(k)=env(ind_xc(it));
                    else
                        tp(k)=nan;
                        ae(k)=nan;
                        ap(k)=nan;
                    end
                end
            case {7}
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    ind_xc=find_zero_crossings(samps2(:,k),0);
                    env=abs(hilbert(samps2(:,k)));
                    if(~isempty(ind_xc))
                        it=near(trel2(ind_xc),0);%find closest to the event time
                        tp(k)=te_int(k)+trel2(ind_xc(it(1)));
                        ae(k)=samps2(ind_xc(it),k);
                        ap(k)=env(ind_xc(it));
                    else
                        tp(k)=nan;
                        ae(k)=nan;
                        ap(k)=nan;
                    end
                end
            case {8}
                if(isempty(THRESH))
                    thresh=2;
                else
                    thresh=THRESH;
                end
                if(isempty(FBAVE))
                    fbave=.05;
                else
                    fbave=FBAVE;
                end
                if(isempty(TIMES))
                    times=10;%the lta will be this many times longer than the sta
                else
                    times=TIMES;
                end
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    env=abs(hilbert(samps2(:,k)));
                    %                     sif=ins_freq(samps2(:,k),trel2);
                    %                     test=env.*sif;
                    %                     bkgrd=mean(test);
                    %                     ind=find(test>thresh*bkgrd);
                    sabs=abs(samps2(:,k));
                    %compute sta/lta
                    n1=max([round(length(env)/(10*times)) 2]);%at least two samples in sta
                    n2=min([times*n1 length(env)]);
                    shortones=ones(1,n1);
                    longones=ones(1,n2);
                    sta=convz(env,shortones);
                    lta=convz(env,longones);
                    ratio=(n2/n1)*sta./lta;
                    if(k==25)
                        honk=1;
                    end
                    %find the average before the event time
                    it=trel2<0;
                    m_early=mean(sabs(it));
                    %we build a purely causal operator of length fbave to
                    %detect the average power of a signal onset.
                    nave=round(fbave/dt2)+1;
                    sabs_ave=convm(sabs,ones(1,nave))/nave;
                    ind=find(sabs_ave>thresh*m_early);
                    if(~isempty(ind)&&ind(1)~=1)
                        %it=near(trel2(ind(1)),0);%find closest to the event time
                        ap(k)=samps2(ind(1),k);
                        ae(k)=samps2(ind(1),k);
                        tp(k)=te_int(k)+trel2(ind(1));
                    else
                        ind=near(trel2,0);
                        %                         ap(k)=env(ind(1));
                        %                         tp(k)=te_int(k)+trel2(ind(1));
                        %                         ae(k)=samps2(ind(1),k);
                        ap(k)=nan;
                        tp(k)=nan;
                        ae(k)=nan;
                    end
                end
            case {9}
                if(isempty(THRESH))
                    thresh=2;
                else
                    thresh=THRESH;
                end
                if(isempty(FBAVE))
                    fbave=.05;
                else
                    fbave=FBAVE;
                end
                %                 if(isempty(TIMES))
                %                     times=10;%the lta will be this many times longer than the sta
                %                 else
                %                     times=TIMES;
                %                 end
                ap=zeros(nx2,1);
                ae=ap;
                tp=ap;
                for k=1:nx2
                    env=abs(hilbert(samps2(:,k)));
%                     sif=ins_freq(samps2(:,k),trel2);
%                     test=env.*sif;
%                     bkgrd=mean(test);
%                     ind=find(test>thresh*bkgrd);
                    [em,ind]=max(env);
                    
                    if(~isempty(ind)&&ind(1)~=1)
                        %it=near(trel2(ind(1)),0);%find closest to the event time
                        ap(k)=samps2(ind(1),k);
                        ae(k)=samps2(ind(1),k);
                        tp(k)=te_int(k)+trel2(ind(1));
                    else
                        ind=near(trel2,0);
                        %                         ap(k)=env(ind(1));
                        %                         tp(k)=te_int(k)+trel2(ind(1));
                        %                         ae(k)=samps2(ind(1),k);
                        ap(k)=nan;
                        tp(k)=nan;
                        ae(k)=nan;
                    end
                end
            otherwise
                error('unknown picking method')
        end
end