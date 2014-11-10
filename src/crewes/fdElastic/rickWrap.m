function [wave,tw] = rickWrap(Dt,fDom)
%[wave,tw] = rickWrap(Dt,fDom)
%Fix up a ricker wavelet for finite-difference use 
% Dt    = The time sample rate in seconds
% fDom  = The central (dominant) frequency
% wave  = The wavelet amplitudes in sequence
% tw    = The wavelet time series in sequence
tlen = Dt*500;
[wavelet,tw1] = ricker(Dt,fDom,tlen);
%plot(wavelet)
%disp(wavelet(1:9))
%stop
[rMax,iMax] = max(wavelet); sc = 1/rMax;
wave1 = sc*wavelet;
%iGood = find(abs(wave1)>.001);
iGood = find(abs(wave1)>.0002);     %Determines taper at end of filter
 %disp(iGood(1))
% figure
% plot(tw1,wave1)
iEnd = iMax*2-iGood(1);
if iEnd>length(wave1)
    error('Time sample rate and wavelet frequency are incompatible')
end
% wave = wave1(iGood(1):len+1-iGood(1));
wave = wave1(iGood(1):iEnd);
len = length(wave);
tw = tw1(1:len);
% disp(rMax); disp(wave(1:5))

