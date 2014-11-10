function [spectrum,f]=multitaper(s,t)
% multitaper: estimate the spectrum of a short seismic trace using multitaper method
%
% [spectrum,f]=multitaper(s,t)
%
% s ... seismic trace
% t ... time coordinate for s
% spectrum ... amplitude spectral estimate
% f ... frequency coordinate for spectrum
%
%
% by Peng Cheng, 2012-2013
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

if(length(s)~=length(t))
    error('s and t must be the same length')
end

nt=length(t);

[sp,f]=fftrl(s,t);
s_est=sp;
k=5;

[e,v]=dpss(nt,4);

spec=zeros(length(sp),2);

for m=1:5
    [spec(:,m),f]=fftrl((s.*e(:,m)),t);
end

spec_est=sum(abs(spec).^2,2)/5;
weight=zeros(size(spec_est));
var=sum(s.^2);
L=20;
for m=1:L
    for n=1:5
        weight(:,n)=(sqrt(v(n))*spec_est)./(v(n)*spec_est+var*(1-v(n)));
    end
    spec_est=sum(abs(spec.*weight).^2,2)./sum(weight,2);
end
tf=isnan(spec_est);
for n=1:length(tf)
    if(tf(n))
        spec_est(n)=0;
    end   
end
spectrum=sqrt(spec_est);

spectrum=spectrum*norm(sp)/norm(spectrum);
