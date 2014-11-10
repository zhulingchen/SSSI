function [wave,tw] = fdMinphase(Dt,fDom,cycleSrc)
%[wave,tw] = fdMinphase(Dt,fDom,cycleSrc)
%Use a double deconned ricker wavelet as a minimum phase wavelet
%The input parameters are
    % Dt    = The time sample rate in seconds
    % fDom  = The central (dominant) frequency
%The output parameters are
    % wave  = The wavelet amplitudes in sequence
    % tw    = The wavelet time series in sequence
tlen = Dt*500;
[wavelet,tw1] = ricker(Dt,fDom,tlen);
len = length(wavelet);
tw2 = tw1-tw1(1);        %Start time series at zero
%Set up shape of minimum phase filter with a double decon
[trout,deconOp] = deconw(wavelet,wavelet,len,.0000001);
[trout,wavelet] = deconw(deconOp,deconOp,len,0);
%Set up taper for wavelet end
A = 1/sqrt(2*pi);
fact = fDom*4.5/(1.25*cycleSrc);
power = ((fact*tw2).^2)*0.5;
env = A*exp(-power);
%Apply taper
wave2 = wavelet.*env;
len = length(wavelet);
%[rMax,iMax] = max(wavelet); sc = 1/rMax;
rMax = max(wave2); sc = 1/rMax;
wave1 = sc*wave2;
%iGood = find(abs(wave1)>.001);
significant = .0002; %.0001;
% iGood = find(abs(wave1)>significant);     %Determines taper at start of filter
%disp(size(wave1))
wave2 = flipud(wave1);
iL = find(abs(wave2)>significant);     %Determines cut off at end of filter
%disp(iL(1:10))
iEnd = len+1-iL(1);
 %disp([len,iGood(1),iEnd])
% % figure
% % plot(tw1,wave1)
% iEnd = iMax*2-iGood(1);
% if iEnd>length(wave1)
%     error('Time sample rate and wavelet frequency are incompatible')
% end
% % wave = wave1(iGood(1):len+1-iGood(1));
%wave = wave1(iGood(1):iEnd);
wave = wave1(1:iEnd);
disp(wave1(iEnd-5:iEnd))
len = length(wave);
tw = tw2(1:len);
% disp(rMax); disp(wave(1:5))

