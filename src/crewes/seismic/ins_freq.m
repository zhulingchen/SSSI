function ifreq = ins_freq(s,t,n,method,flag)
% INS_FREQ Compute the instantaneous frequency using complex trace theory
%
% ifreq=ins_freq(s,t,n,method,flag)
%
% See Barnes, 1992, Geophysics, vol 57, November, P. 1520-1524
%
% s ... real seismic trace (or gather)
% t ... time coordinate vector for s (or time sample interval)
% n ... length of a boxcar smoother (in samples) to be convolved with the result
%  ************* default n=1 **********
% method ... 1 Unwrap the instantaneous phase and differentiate (Barnes eqn 2)
%            2 Differentiate real and imaginary parts (Barnes Eqn 3)
%            3 Claerbout's method (Barnes eqn 14)
%            4 Barnes first approximation (Barnes eqn 9)
%            5 Barnes second approximation (Barnes eqn 12)
% ******* default method=2 ********
% flag ... 1 Calculate the absolute value of the instantaneous frequency
%      ... 2 Calculate the signed instantaneous frequency
% *********** default = 1 *********
%
% ifreq ... instantaneous frequency of the complex trace corresponding to s
% NOTE s may be a single trace or a gather.
%
% G.F. Margrave, CREWES, 2014
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

if(nargin<3)
    n=1; 
end
if(nargin<4)
    method=2;
end
if(nargin<5)
    flag=1;
end
if(size(s,1)==1)
    s=s';%avoid row vectors
end


if(method==1)
    %classic way
    iph=ins_phase(s);
    ifreq=zeros(size(iph));
    ntr=size(iph,2);
    if(length(t)>1)
        dt=t(2)-t(1);
    else
        dt=t;
    end
    oneupontwopi=1/(2*pi);
    for k=1:ntr
        ifreq(:,k)=oneupontwopi*gradient(unwrap(iph(:,k)),dt);
    end
elseif(method==2)
    %Barnes exact
    ifreq=zeros(size(s));
    ntr=size(s,2);
    if(length(t)>1)
        dt=t(2)-t(1);
    else
        dt=t;
    end
    oneupontwopi=1/(2*pi);
    for k=1:ntr
        sa=hilbert(s(:,k));%analytic trace
        x=real(sa);
        y=imag(sa);
        xp=gradient(x,dt);
        yp=gradient(y,dt);
        ifreq(:,k)=oneupontwopi*(x.*yp-xp.*y)./(x.^2+y.^2);
    end
elseif(method==3)
    %claerbouts method
    ifreq=zeros(size(s));
    ntr=size(s,2);
    if(length(t)>1)
        dt=t(2)-t(1);
    else
        dt=t;
    end
    twouponpidt=2/(pi*dt);
    for k=1:ntr
        sa=hilbert(s(:,k));%analytic trace
        x=real(sa);
        y=imag(sa);
        x1=[x(2:end);0];
        y1=[y(2:end);0];
        ifreq(:,k)=twouponpidt*(x.*y1-x1.*y)./((x+x1).^2+(y+y1).^2);
    end
elseif(method==4)
    %Barnes first approximation
    ifreq=zeros(size(s));
    ntr=size(s,2);
    if(length(t)>1)
        dt=t(2)-t(1);
    else
        dt=t;
    end
    oneupontwopidt=1/(2*pi*dt);
    for k=1:ntr
        sa=hilbert(s(:,k));%analytic trace
        x=real(sa);
        y=imag(sa);
        x1=[x(2:end);0];
        y1=[y(2:end);0];
        ifreq(:,k)=oneupontwopidt*atan2(x.*y1-x1.*y,x.*x1+y.*y1);
    end
elseif(method==5)
    %Barnes second approximation
    ifreq=zeros(size(s));
    ntr=size(s,2);
    if(length(t)>1)
        dt=t(2)-t(1);
    else
        dt=t;
    end
    oneuponfourpidt=1/(4*pi*dt);
    for k=1:ntr
        sa=hilbert(s(:,k));%analytic trace
        x=real(sa);
        y=imag(sa);
        x1=[x(2:end);0];
        y1=[y(2:end);0];
        x2=[0;x(1:end-1)];
        y2=[0;y(1:end-1)];
        ifreq(:,k)=oneuponfourpidt*atan2(x2.*y1-x1.*y2,x2.*x1+y2.*y1);
    end
end
if(flag==1)
    ifreq=abs(ifreq);
end
if(n>1)
    ifreq=convz(ifreq,ones(n,1))/n;
end
    
