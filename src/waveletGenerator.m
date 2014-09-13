% Matlab script for the computation of the fdmpi source signal
%
% Daniel Koehn
% Freiberg, 6th of August 2008

close all;
clear;
clc;

% choose wavelet type
% 1 = Ricker wavelet (as in fdmpi)
% 2 = Fuchs-Mueller wavelet
% 3 = SINE^3 wavelet
% 4 = Ricker wavelet (from H. Ryan, "Ricker, Ormsby, Klauder, Butterworth - A Choice of Wavelets", CSEG Recorder, Sep. 1994)
% 5 = Ormsby wavelet
% 6 = Klauder wavelet

for wavelet = 1:6
    
    % Define Parameters from fdmpi input file
    dt = 5.4e-4;    % timestep [s]
    T1 = 0.0;      % first temporal sample point [s]
    T2 = 1.5;     % time of wave propagation [s]
    fc = 10.0;    % center frequency of the source wavelet [Hz]
    tshift = 0.8; % time shift of the source signal
    amp = 1.0;    % amplitude of the source signal
    
    % some internal parameters
    tsour=1/fc;      % source length [s]
    t=[T1:dt:T2];    % define time
    n=length(t);     % number of samples
    lsour=tsour/dt;  % length of the source wavelet [sample points]
    
    % plotting parameters for the amplitude spectrum
    % plot frequencies between f1 and f2 Hz
    fp1=-140.0;
    fp2=140.0;
    
    % parameters for the Ormsby wavelet
    fo1 = 5;
    fo2 = 10;
    fo3 = 40;
    fo4 = 45;
    
    % parameters for the Klauder wavelet
    fk1 = 10;
    fk2 = 40;
    Tk = 6.0;
    
    % calculate source signal for fdmpi
    
    if (wavelet==1)
        % RICKER-SIGNAL:
        t0=tsour;
        T0=tsour*1.5;
        tau=pi.*(t-t0-tshift)/T0;
        ft=(1.0-4.0*tau.*tau).*exp(-2*tau.*tau);
    end
    
    if (wavelet==2)
        % FUCHS-MUELLER-SIGNAL
        w1=2*pi.*fc;
        
        for ii=1:n
            
            ft(ii)=sin(w1.*(t(ii)-tshift))-0.5*sin(2*w1*(t(ii)-tshift));
            
            if((t(ii)<tshift)||(t(ii)>(tshift+tsour)))
                ft(ii) = 0.0;
            end
            
        end
    end
    
    if (wavelet==3)
        % SINE^3-SIGNAL
        
        for ii=1:n
            
            ft(ii)=(0.75.*pi.*fc).*(sin(pi.*(t(ii)+tshift).*fc)).^3.0;
            if((t(ii)<tshift)||(t(ii)>(tshift+tsour)))
                ft(ii) = 0.0;
            end
            
        end
    end
    
    if (wavelet==4)
        % RICKER-SIGNAL (modified):
        tau=pi.*(t-tshift).*fc;
        ft=(1.0-2.0*tau.*tau).*exp(-2*tau.*tau);
    end
    
    if (wavelet==5)
        % ORMSBY-SIGNAL:
        tau1 = ((pi.*fo4).^2)./((pi.*fo4)-pi.*fo3);
        tau2 = ((pi.*fo3).^2)./((pi.*fo4)-pi.*fo3);
        tau3 = ((pi.*fo2).^2)./((pi.*fo2)-pi.*fo1);
        tau4 = ((pi.*fo1).^2)./((pi.*fo2)-pi.*fo1);
        
        % find sample point at which t-tshift is a minimum
        ts = abs(t-tshift);
        [minshift,ishift] = min(ts);
        
        % calculate sinc function
        x = fo4.*(t-tshift);
        sinc1 = sin(pi.*x)./(pi.*x);
        sinc1(ishift) = 1.0;
        
        x = fo3.*(t-tshift);
        sinc2 = sin(pi.*x)./(pi.*x);
        sinc2(ishift) = 1.0;
        
        x = fo2.*(t-tshift);
        sinc3 = sin(pi.*x)./(pi.*x);
        sinc3(ishift) = 1.0;
        
        x = fo1.*(t-tshift);
        sinc4 = sin(pi.*x)./(pi.*x);
        sinc4(ishift) = 1.0;
        
        % calculate Ormsby signal
        ft = ((tau1 .* sinc1.^2) - (tau2 .* sinc2.^2)) - ((tau3 .* sinc3.^2) - (tau4 .* sinc4.^2))
        
        ft = ft ./max(ft);
        
    end
    
    if (wavelet==6)
        % KLAUDER-SIGNAL:
        k = (fk2-fk1)./Tk;
        f0 = (fk2+fk1)./2.0;
        ft = real((sin(pi.*k.*(t-tshift).*(Tk-(t-tshift)))./(pi.*k.*(t-tshift))).*exp(2.*pi.*1j.*f0.*(t-tshift)));
        ft = ft ./max(ft);
    end
    
    ft=amp.*ft;
    
    % Plot of the source signal in the time domain
    figure;
    plot(t,ft)
    title(['Source Signal']);
    xlabel('Time [s]');
    ylabel('Amplitude');
    
    fnyq=1/(2*dt);
    nfft=2^nextpow2(length(ft));
    
    %FFT
    y=fftshift(fft(ft,nfft));
    amp=abs(y);
    amp=amp/max(amp);
    f=fnyq*(-nfft/2:nfft/2-1)/(nfft/2);
    
    nn1=round((nfft/2+1)+(fp1/fnyq)*nfft/2);
    nn2=round((nfft/2+1)+(fp2/fnyq)*nfft/2);
    
    % Plot of the source signal in the frequency domain
    figure;
    plot(f(nn1:nn2),amp(nn1:nn2))
    title('Amplitude-Spectrum')
    xlabel('Frequency [Hz]')
    ylabel('Amplitude');
    
end