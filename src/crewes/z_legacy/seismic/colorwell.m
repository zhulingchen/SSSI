function [cop,tcop,specmodel,fsm]=colorwell(r,t,n)
% COLORWELL: compute an operator to restore well color after decon
%
% [cop,tcop,specmodel,fsm]=colorwell(r,t,n)
%
% COLORWELL designs a convolutional minimum phase operator which
% can be applied to the output of any spiking decon to correct
% for a non-white reflectivity.
%
% r= well log reflectivity series whose amplitude spectrum shows the expected 
%    non-white behavior
% t= time coordinate for r
% n= desired length of the color operator (samples)
%
% cop= output color operator. This will be minimum phase and of length n
% tcop= time coordinate for cop
% specmodel= spectral model designed to approximate reflectivity spectrum
%       at low frequencies
% fsm= frequency coordinate for specmodel
%
% Compare the specmodel to the well with
% [R,f]=fftrl(r,t)
% plot(f,abs(R),fsm,specmodel)
% 
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
  
%compute reflectivity spectrum
r=padpow2(r);
t=(t(2)-t(1))*(0:length(r)-1);
[R,f]= fftrl(r,t);
% fit a curve of the form a+b*atan(f/f0) to the spectrum from 1Hz to 1/2 Nyquist
% Here f0 will be taken as 1/4 Nyquist
fhalfnyq=f(end)/2;
% f0=f(end)/10;
%search for optimal f0
ftest=40:120;
obj=zeros(size(ftest));
indf=near(f,1,fhalfnyq);
for k=1:length(ftest)
    f0=ftest(k);
    A=[ones(size(f(indf))) atan(f(indf)/f0)];
    ab=A\abs(R(indf));%solve the spectral model
    specmodel=zeros(size(f));
    specmodel(2:end)=ab(1)+ab(2)*atan(f(2:end)/f0);
    specmodel(1)=specmodel(2);%avoid a zero at zero Hz
    obj(k)=sum((specmodel(1:ftest(end))-abs(R(1:ftest(end)))).^2);
end
[omin,imin]=min(obj);
f0=ftest(imin);
A=[ones(size(f(indf))) atan(f(indf)/f0)];
ab=A\abs(R(indf));%solve the spectral model
specmodel=zeros(size(f));
specmodel(2:end)=ab(1)+ab(2)*atan(f(2:end)/f0);
specmodel(1)=specmodel(2);%avoid a zero at zero Hz
% compute the minimum phase time domain operator
cop= minwave(specmodel,f,0,n);
fsm=f;
tcop=(0:length(cop)-1)*(t(2)-t(1));
