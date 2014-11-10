function qmat=qmatrix(q,t,w,tw,flag)
%QMATRIX ... create a matrix that applies a forward Q filter
%
% qmat=qmatrix(q,t,w,tw,flag)
%
% Generates a matrix which applies a forward Q filter. The qmatrix is also
% bandlimited by a stationary "source signature" waveform. Each column of
% the Q matrix is the theoretical Q impulse response for the traveltime of
% the column convolved with the stationary source waveform. The theoretical
% Q impulse response is computed by einar.m and is based on Kjartanssen
% (1979). The group delay induced by attenuation is often called the drift.
% This is controlled in einar.m by the sample-rate at which the impulse
% response is calculated.
%
% q ... value of Q, may be a single scalar or a vector of length(t)
%        of average Q values
% t ... vector of times to generate q response at
% w ... stationary source waveform
% tw ... time vector for w
% *********** Caution ************
% Q matrix only works properly with causal wavelets. Ideally w should be
% minimum phase.
% *********** Caution ************
% qmat ... length(t) by length(t) matrix which applies a forward Q
%		filter.
% 		If r is a column vector of reflection coefficients, then
%		s=qmat*r; is a simple synthetic seismogram with attenuation.
% flag ... 1 means calcuate the Q impulse response at the seismic sample
%               rate (tw(2)-tw(1))
%      ... 2 means calcuate the Q impulse response at a sample rate of
%           .00002 s whach has a Nyquist of 25kHz and the resample to the
%           seismic sample rate. This simulates the group delay relative to
%           a logging frequency of 12.5kHz but is very expensive.
%      ... 3 means calculate the Q impulse response at the seismic sample
%            rate but then delay it using the theoretical drift time as
%            calculated by tdrift.m .
%     Remark: flag=2 takes a long time. Flag=3 is much quicker and produces
%     similar results. Flag=1 is the simplest option but get unrealistically
%     small drift times.
%     ************** default flag=1 ****************
%
% Example (Copy and paste these lines to the Matlab command window)
% dt=.002;%time sample rate
% tmax=2;%maximum trace time
% Q=50;%Q value
% fdom=30;%dominant frequency of minimum phase source waveform
% [r,t]=reflec(tmax,dt);%synthetic reflectivity
% [w,tw]=wavemin(dt,fdom,.2);
% qmat=qmatrix(Q,t,w,tw);%build the qmatrix
% sn=qmat*r;%apply the Q matrix to the reflectivity
% s=convm(r,w);%stationary synthetic for comparison
% figure
% plot(t,s,t,sn,'r')
% xlabel('time (s)')
% legend('stationary',['nonstationary, Q=' num2str(Q)])
%
% G.F. Margrave, July 1996,2013 The CREWES Project
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
    flag=1;
end
if(flag<1 || flag>3)
    error('invalid value for flag')
end

v=2000;%a fake velocity that is not relevant in the end.
x=v*t;
nt=length(t);
tmax=max(t);
dt=t(2)-t(1);
fnyq=1/(2*dt);
dt2=.00002;
fnyq2=1/(2*dt2);
qmat=zeros(nt,nt);
if(nargin>2)
	izero=near(tw,0);
    fd=dom_freq(w,tw);
end
if(length(q)==1)
    q=q*ones(size(t));
elseif(length(q)~=length(t))
    error('q must be a scalar or a vector of length(t), no other options')
end
if(flag~=1)
    td1=tdrift(q,x,fd,v*ones(size(x)),fnyq/2);
    td2=tdrift(q,x,fd,v*ones(size(x)),fnyq2/2);
end

for k=1:nt
    if(flag==1)
        tmp=einar(q(k),x(k),v,dt,tmax,1);
    elseif(flag==2)
        [tmp2,t2]=einar(q(k),x(k),v,dt2,tmax,1);
        tmp=resamp(tmp2,t2,dt,[0 tmax],1);
    else
        [tmp2,t2]=einar(q(k),x(k),v,dt,tmax,1);
        %now delay by the difference between the drift times at fnyq/2 and
        %fnqy2/2
        tmp=stat(tmp2,t2,td2(k)-td1(k));
    end
% 	nuse=round((tmax-t(k))/dt) +1;
	
	%include stationary wavelet
	if(nargin>2)
		tmp=convz(tmp,w,izero,nt,0);
	end
% 	if(k==1)
%         amax=max(abs(tmp));
%     end
% 	qmat(k:nt,k)=tmp(1:nuse)'/amax;
%     qmat(:,k)=tmp'/amax;
    qmat(:,k)=tmp';
end

%wconv=convmtx(w(:),nt);
%qmat=wconv(1:nt,:)*qmat;
	
	
