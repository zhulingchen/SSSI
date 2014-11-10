function Qint=qestimator(s1,s2,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true)
%
% Qint=qestimator(s1,s2,t,t1,t2,twin,f1,f2,wintype,spectype,pflag,method,Qint_true)
%
% s1 ... trace 1 (trace with less attenuation)
% s2 ... trace 2 (trace with greater attenuation)
% NOTE: s1 and s2 can be the same trace but t2 should be greater than t1 so
%   that s2 at t2 shows more attenuation than s1 at t1.
% t  ... time coordinate for traces
% t1 ... start time (scalar) for spectal estimate on trace 1
% t2 ... start time (scalar) for spectal estimate on trace 2
% twin ... temporal window length for spectral estimation
% f1 ... lowest frequency to allow in spectral ratio least-square fit
% f2 ... highest frequency to allow in spectral ratio least-squares fit
% wintype ... 'boxcar' 
%             'cosine' (raised cosine)
%             'gaus' Gaussian truncated at two standard deviations
%             'mwin'  see mwindow.m (20% taper)
%             'mwhalf' see mwhalf.m (20% taper). Intended for VSPs
% NOTE: wintype only matters if spectype='fourier'
% spectype ... 'fourier' Fourier transform of windowed signal
%              'burg' Burg spectrum of boxcar widowed signal segment
%              'multi' Multitaper spectral estimate of boxcar windowed signal segment
% pflag ... plot flag, 1 means plot, 0 means don't plot
% method ... 'specrat' Spectral ratio method
%            'specmatch' Amplitude spectral matching
%            'domfreq' Dominant frequency matching
%            'all' Return the median from all methods
% Qint_true ... this is only used for plot annotation. It has no effect on
%           the estimation. Omit it if you wish.
%
% Qint ... estimated interval Q between t1 and t2
%
% by G.F. Margrave, 2013-2014
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

if(nargin<13)
    Qint_true=[];
end
if(length(s1)~=length(s2))
    error('s1 and s2 must be the same length')
end
if(length(t)~=length(s1))
    error('t and s1 must be the same length');
end
if(t2<t1)
    error('s1 and t1 should refer to the trace with less attenuation, but you have t2<t1')
end

if(~strcmp(wintype,'boxcar')&&~strcmp(wintype,'cosine')&&~strcmp(wintype,'gaus')...
        &&~strcmp(wintype,'mwin')&&~strcmp(wintype,'mwhalf'))
    error('Invalid wintype')
end
if(~strcmp(spectype,'fourier')&&~strcmp(spectype,'burg')&&~strcmp(spectype,'multi'))
    error('Invalid spectype')
end
if(~strcmp(method,'specrat')&&~strcmp(method,'specmatch')&&~strcmp(method,'domfreq')&&~strcmp(method,'all'))
    error('Invalid method')
end


dt=t(2)-t(1);
t1=round(t1/dt)*dt;%make sure it falls on a sample
t2=round(t2/dt)*dt;%ditto
delt=t2-t1;
%build window
ntwin=round(twin/dt)+1;
switch wintype
    case 'boxcar'
        win=ones(ntwin,1);
        winname='boxcar';
    case 'cosine'
        x=(1:ntwin)';
        win=.5+.5*cos(-pi*(x-ntwin)/(1-ntwin)+pi*(x-1)/(ntwin-1));
        winname='raised cosine';
    case 'gaus'
        x=(1:ntwin)';
        x0=(1+ntwin)/2;
        %want a gaussian centered at x0 with sigma of ntwin/4
        sigma=ntwin/4;
        win=exp(-(x-x0).^2/sigma^2);
        winname='gaussian';
    case 'mwin'
        pct=20;
        win=mwindow(ntwin,pct);
        winname='mwindow';
    case 'mwhalf'
        pct=20;
        win=mwhalf(ntwin,pct);
        winname='mwhalf';
        
    otherwise
        error('invalid wintype')
end

nmethods=3;

Qintall=nan*ones(1,nmethods);

%compute spectra
inwin1=near(t,t1,t1+twin);
inwin2=near(t,t2,t2+twin);
if(strcmp(spectype,'fourier'))
    sw1=win.*s1(inwin1);
    sw2=win.*s2(inwin2);
    [S1,f]=fftrl(sw1,t(inwin1));
    S2=fftrl(sw2,t(inwin2));
elseif(strcmp(spectype,'burg'))
    sw1=s1(inwin1);
    sw2=s2(inwin2);
    L=min([round(length(inwin1)/10) 10]);
%     L=20;
    [S1,f]=burg(sw1,t(inwin1),L);
    S2=burg(sw2,t(inwin2),L);
elseif(strcmp(spectype,'multi'))
    sw1=s1(inwin1);
    sw2=s2(inwin2);
    [S1,f]=multitaper(sw1,t(inwin1));
    S2=multitaper(sw2,t(inwin2));
else
    error('invalid spectype')
end

if(strcmp(method,'specrat')||strcmp(method,'all'))
%     %spectral ratio
    fsmo=0;
    [Q,T,SR,SR_fit,f1e,f2e]=q_specrat(abs(S1),abs(S2),f,t1,t2,f1,f2,fsmo);
    infit=near(f,f1,f2);
    SR_fit=SR_fit(infit);
    Qintall(1)=Q;
end
if(strcmp(method,'specmatch')||strcmp(method,'all'))
    %spectral matching
    Qmin=5;
    Qmax=250;
    [Q,T,obj]=q_specmatch(abs(S1),abs(S2),f,t1,t2,f1,f2,Qmax,Qmin);
    Qintall(2)=Q;
    Qtest=Qmin:Qmax;
    iq=near(Qtest,Q);
end
if(strcmp(method,'domfreq')||strcmp(method,'all'))
%     Qtest=5:250;
%     A1=abs(S1);
%     A2=abs(S2);
%     fd1=sum(f.*A1.^2)/sum(A1.^2);
%     fd2=sum(f.*A2.^2)/sum(A2.^2);
%     if(fd1>fd2)
%         fdtest=zeros(size(Qtest));
%         for k=1:length(Qtest)
%             Atest=A1.*exp(-pi*f*delt/Qtest(k));
%             fdtest(k)=sum(f.*Atest.^2)/sum(Atest.^2);
%         end
%         %find where fdtest is closest to fd2;
%         obj=Qtest.*(fdtest-fd2).^2;%objective function
%         [omin,imin]=min(obj);
%         Qintall(3)=Qtest(imin);
%     end
    %centroid frequency matching
    Qmin=5;
    Qmax=250;
    [Q,obj,fd1,fd2]=q_centroid(abs(S1),abs(S2),f,t1,t2,Qmax,Qmin);
    Qintall(3)=Q;
    Qtest=Qmin:Qmax;
    iq=near(Qtest,Q);
end

if(pflag==1&&~strcmp(method,'all'))
    figure
    if(strcmp(method,'specmatch'))
        pos=get(gcf,'position');
        set(gcf,'position',[pos(1) pos(2)/1.25 pos(3) 1.25*pos(4)]);
        npanes=4;
        if(strcmp(spectype,'fourier'))
            titl=['Spectral matching (Fourier) using ' winname ' windows'];
        elseif(strcmp(spectype,'burg'))
            titl='Spectral matching using Burg spectra';
        else
            titl='Spectral matching using Multitaper spectra';
        end
    elseif(strcmp(method,'specrat'))
        npanes=3;
        if(strcmp(spectype,'fourier'))
            titl=['Spectral ratio (Fourier) estimate using ' winname ' windows'];
        elseif(strcmp(spectype,'burg'))
            titl='Spectral ratio using Burg spectra';
        else
            titl='Spectral ratio using Multitaper spectra';
        end
    elseif(strcmp(method,'domfreq'))
        npanes=3;
        if(strcmp(spectype,'fourier'))
            titl=['Dominant frequency (Fourier) estimate using ' winname ' windows'];
        elseif(strcmp(spectype,'burg'))
            titl='Dominant frequency using Burg spectra';
        else
            titl='Dominant frequency using Multitaper spectra';
        end
    end
    subplot(npanes,1,1)
    shift=max(abs(sw1));
    plot(t(inwin1),sw1,t(inwin2),sw2+shift,'r')
    title(titl)
    xlabel('time (s)')
    legend(['s1 (twin=' num2str(twin) 's)'],'s2')
    subplot(npanes,1,2)

    if(strcmp(method,'domfreq'))
        plot(f,abs(S1),f,abs(S2),'r')
        xlabel('frequency (Hz)')
        legend(['S1 (fdom= ' num2str(fd1) 'Hz) '],...
        ['S2 (fdom= ' num2str(fd2) 'Hz) ']);
    else
        plot(f,abs(S1),f,abs(S2),'r',f,abs(S2)/T)
        xlabel('frequency (Hz)')
        legend(['S1 (t1= ' num2str(t1) 's) '],...
        ['S2 (t2= ' num2str(t2) 's) '],['S2/T (T= ' num2str(T) ')']);
    end
    
    if(strcmp(method,'specrat'))
        subplot(npanes,1,3)
        h1=plot(f,SR);
        xlabel('frequency (Hz)')
        hq=line(f(infit),SR_fit,'color','r');
        ymin=min(SR_fit)-abs(min(SR_fit));
        ymax=max(SR_fit)+abs(max(SR_fit));
        xl=get(gca,'xlim');
        ylim([ymin ymax])
        xlim(xl);
        if(isempty(Qint_true))
             legend([h1 hq],'LSR= log(S_2/S_1)',['Best fit line, Qint=' int2str(Qintall(1))])
        else
            legend([h1 hq],'LSR= log(S_2/S_1)',['Best fit line, Qint=' int2str(Qintall(1))...
            ' actual =' int2str(Qint_true)])
        end
    elseif(strcmp(method,'specmatch'))
        subplot(npanes,1,3)
        plot(f,T*abs(S1).*exp(-pi*f*delt/Qintall(2)),f,abs(S2))
        xlabel('frequency (Hz)')
        if(isempty(Qint_true))
             legend('R*S_1.*exp(-pi*f*delt/Qint)',['S_2, Qint=' int2str(Qintall(2))])
        else
            legend('R*S_1.*exp(-pi*f*delt/Qint)',['S_2, Qint=' int2str(Qintall(2)) ...
            ' actual =' int2str(Qint_true)])
        end
        subplot(npanes,1,4)
        err_min=nan*obj;
        err_min(iq)=obj(iq);
        plot(Qtest,obj,Qtest,err_min,'r*')
        xlabel('Q value');ylabel('Error norm')
        legend('Error norm','Optimal value')
    elseif(strcmp(method,'domfreq'))
        subplot(npanes,1,3)
        ominn=nan*obj;
        ominn(iq)=obj(iq);
        plot(Qtest,obj,Qtest,ominn,'r*')
        xlabel('Q value')
        legend('objective function',['minimimum, Q_{int}=' num2str(Qintall(3))])
    end
    
end

if(strcmp(method,'all'))
    ind=find(~isnan(Qintall));
    if(~isempty(ind))
        Qint=median(Qintall(ind));
    else
        Qint=nan;
    end
elseif(strcmp(method,'specrat'))
    Qint=Qintall(1);
elseif(strcmp(method,'specmatch'))
    Qint=Qintall(2);
elseif(strcmp(method,'domfreq'))
    Qint=Qintall(3);
end