%define the synthetic to be deconvolved
%A random reflectivity is generated with "reflec". The frist time you run
%this script in a session, the synthetic attenuated trace is built. All
%subsequent times, the existing trace is used. This allows the Gabor
%parameters to be studied against the same underlying reflectivity.  The
%script checks to see if the variable rq exists, which is the random
%reflectivity. If so, then a new synthetic is not made regardless of what
%you may have changed q or tmax to. To force a new synthetic to be made,
%you must clear the variable rq: >>clear rq
%If you wish to terminate your Matlab session and restart later with the
%same synthetic, just use the save command to save your workspace and then
%load it in the new Matlab session before running this script
dt=.002;tmax=2;q=50;

iburg=0;%set to zero for Fourier algorithm... Burg not currently recommended
iwiener=1;%set to 1 for an AGC-WIENER-AGC comparison, leave this alone or risk trouble
iopt=1;%1 means analysis windows, 2 is synthesis, 3 is both
ibigfigs=0;% set to zero to turn of big figures

idemo=0;%demo forward and inverse gabor
idecon=1;%deconvolve or not
ibandlimit=0;% 0  means don't filter after decon, 0 means filter
igaborfact=0;% set to zero to turn off plotting of the factors of the input Gabor spectrum

%gabordecon parameters, type >>help gabordecon 
%for more information
twin=.2;tinc=.01;%defines the Gaussian windows
gdb=60;%truncation factor in fgabor
pow2option=1;%make windows a power of 2 in length
%the next three parameters define the smoothing operation on the Gabor
%spectrum of the input signal
tsmo=.4;%temporal smoother in seconds
fsmo=2;%frequency smoother in Hz
ihyp=1;%flag for hyperbolic smoothing. 1 gets Hyperbolic, 0 gets boxcar
%Define the p parameter that controls the tradeoff between analysis and
%synthesis windows. p=1 gets all analysis (preferred) p=0 gets all
%synthesis (bad)
p=1;
%
stabg=0;phase=1;%gabor decon stab factor and phase flag
order=10;%order of the Burg spectrum if iburg is 1

%wiener decon parameters
stabw=.0001;%stability constant
operator_length=100;%operator length in samples



%plotflags
plottraces=1; %plot time domain traces
plotfourier=0; %plot Fourier amplitude spectra of traces
gaborplot=1;% set to zero to turn off all ploting of Gabor spectra regardless of the next flags
plotgab_input=0; %plot Gabor spectrum of input (attenuated) signal
plotgab_wavelet=0; %plot Gabor spectrum of estimated wavelet
plotgab_refl=0; %plot Gabor spectrum of actual reflectivity
plotgab_refl_est=1; %plot Gabor spectrum of reflectivity estimate
plotgab_wiener=0;%plot Gabor spectrum of Wiener estimate of reflectivity

if(iburg)
    titalg=' Burg,';
else
    titalg=' Fourier,';
end
if(iopt==1)
    titwin=' analysis,';
elseif(iopt==2)
    titwin=' synthesis,';
elseif(iopt==3)
    titwin=' analysis and synthesis,';
end
if(~exist('rq'))
    [rq,t]=reflec(tmax,dt,.2);
    [w,tw]=wavemin(dt,20,tmax/10);
    qmat=qmatrix(q,t,w,tw);
    %s=convm(r,w);
    sq=qmat*rq;
end

if(idemo) %demo forward and inverse transforms
    normflag=0;
   [sqtvs,trow,fcol]=fgabor(sq,t,twin,tinc,p,gdb,normflag);

   sqq=igabor(sqtvs,trow,fcol,twin,tinc,p,gdb,normflag);
   
   sqq=sqq(1:length(sq));%unpad
   
   figure;
   inc=.8*max(sq);

   subplot(2,1,1)
   plot(t,sq,t,sqq+inc);
   xlabel('Time in seconds')
   title('signal reconstruction after forward and inverse gabor')
   legend('Attenuated signal','After forward and inverse Gabor');
   subplot(2,1,2)
   plot(t,sq-sqq);
   title('difference')
%     nudge=-.1*inc;
%     start=.5*tmax;
%     text(start,nudge,'Attenuated signal');
%     text(start,nudge+inc,'After forward+inverse Gabor');
%     text(start,nudge+2*inc,'Difference');
%     yoff
    
    if(ibigfigs)
        bigfig;whitefig;boldlines(gca,6);bigfont;
    end
end

if(idecon)
    if(~iburg)
       [r2,tvs_op]=gabordecon(sq,t,twin,tinc,tsmo,fsmo,ihyp,stabg,phase,p,gdb,transforms);
    else
       [r2,tvs_op]=gabordeconb(sq,t,twin,tinc,tsmo,fsmo,ihyp,order,stabg,phase,p,gdb);
       %r2=gabordeconbq(sq,t,twin,tinc,tsmo,fsmo,10,stab,phase,q,iopt);
    end


    if(iwiener)
        % wiener for comparison
        izone=near(t,.3*tmax,.7*tmax);
        sqa=aec(sq,t,tmax/3);
        r3=deconw(sqa,sqa(izone),operator_length,stabw);
        r3=aec(r3,t,tmax/3);
        r3=balans(r3,rq);
    end

    fnyq=1/(2*dt);
    fmax=fnyq/2;
    fwid=fnyq/20;
    if(plottraces)
        if(ibandlimit==1)
            r2f=filtf(r2,t,[0 0],[fmax fwid]);
            r2f=balans(r2f,rq);
            rqf=filtf(rq,t,[0 0],[fmax fwid]);
            rqf=balans(rqf,rq);
            r3f=filtf(r3,t,[0 0],[fmax fwid]);
            r3f=balans(r3f,rq);
        else
            r2f=balans(r2,rq);
            rqf=rq;
            r3f=balans(r3,rq);
        end
        figure;
        %subplot(2,1,1)
        sq2=balans(sq,rq);
        plot(t,sq2,'b',t,r3f+.1,'g',t,r2f+.2,'r',t,rqf+.3,'k')
        xlabel('Time in seconds')
        title([titalg titwin 'twin=' num2str(twin) ...
                ' tinc=' num2str(tinc) ' tsmo=' num2str(tsmo) ...
                ' fsmo=' num2str(fsmo) ' ihyp=' num2str(ihyp) ...
                ' stab=' num2str(stabg)])
        nudge=.05;
        start=.7*tmax;
        text(start,nudge,'Attenuated signal');
        text(start,nudge+.1,'After AGC+Wiener');
        text(start,nudge+.2,['After Gabor ' titalg(1:end-1)]);
        text(start,nudge+.3,'True Reflectivity');
        yoff
        if(ibandlimit)
            title(['Estimates bandlimited to ' num2str(fmax) ' Hz']);
        else
            title(['Broadband estimates']);
        end

        if(ibigfigs)
                bigfig;whitefig;boldlines(gca,6);bigfont;
        end
    end
    %subplot(2,1,2)
    if(plotfourier)
        figure
        [R,f]=fftrl(rq2,t);
        R2=fftrl(r2,t);
        R3=fftrl(r3,t);
        SQ=fftrl(sq2,t);
        plot(f,abs(SQ)/max(abs(SQ)),'b',f,abs(R3)/max(abs(R3))+1,'g',f,abs(R2)/max(abs(R2))+2,'r',f,abs(R)/max(abs(R))+3,'k')
        xlabel('Frequency in Hz')
        title([titalg titwin 'twin=' num2str(twin) ...
                ' tinc=' num2str(tinc) ' tsmo=' num2str(tsmo) ...
                ' fsmo=' num2str(fsmo) ' stab=' num2str(stab)])
        nudge=.2;
        start=.6*fnyq;
        text(start,nudge,'Attenuated signal');
        text(start,nudge+1,'After AGC+Wiener+AGC');
        text(start,nudge+2,['After Gabor ' titalg(1:end-1)]);
        text(start,nudge+3,'Bandlimited reflectivity');
        yoff

        if(ibigfigs)
                bigfig;whitefig;boldlines(gca,6);bigfont;
        end
    end

    %make some gabor plots
    if(gaborplot==1)
       if(plotgab_input)

           [sqtvs,trow,fcol]=fgabor(sq,t,twin,tinc,1,gdb);
           plotimage(abs(sqtvs),trow,fcol);
           title('Attenutated signal');
           xlabel('Frequency in Hz');ylabel('Time in seconds')
           if(ibigfigs)
            bigfig;whitefig;boldlines;bigfont;
           end
       end
       if(plotgab_wavelet)
           if(~exist('trow'))
               [sqtvs,trow,fcol]=fgabor(sq,t,twin,tinc,1,gdb);
           end
           plotimage(1./abs(tvs_op),trow,fcol);
           title('Propagating wavelet');
           xlabel('Frequency in Hz');ylabel('Time in seconds')
           if(ibigfigs)
            bigfig;whitefig;boldlines;bigfont;
           end
       end

       if(plotgab_refl)

           [rqtvs,trow,fcol]=fgabor(rq,t,twin,tinc,1,gdb);
           plotimage(abs(rqtvs),trow,fcol);
           title('Reflectivity');
           xlabel('Frequency in Hz');ylabel('Time in seconds')
           if(ibigfigs)
            bigfig;whitefig;boldlines;bigfont;
           end
       end
       if(plotgab_refl_est)

           [r2tvs,trow,fcol]=fgabor(r2,t,twin,tinc,1,gdb,pow2option);
           plotimage(abs(r2tvs),trow,fcol);
           if(iburg)
            title('After Gabor Decon (Burg)');
           else
            title('After Gabor Decon (Fourier)');
           end
           xlabel('Frequency in Hz');ylabel('Time in seconds')
           if(ibigfigs)
            bigfig;whitefig;boldlines;bigfont;
           end
       end

       if(plotgab_wiener)

           [r3tvs,trow,fcol]=fgabor(r3,t,twin,tinc,1,gdb);
           plotimage(abs(r3tvs),trow,fcol);
           title('After AGC+Wiener Decon');
           xlabel('Frequency in Hz');ylabel('Time in seconds')
           if(ibigfigs)
            bigfig;whitefig;boldlines;bigfont;
           end
       end

     end

     if(igaborfact)
       alpha=exp(-pi*trow(:)*fcol/q);
       plotimage(alpha,trow,fcol);
       title('Constant Q surface');
       xlabel('Frequency in Hz');ylabel('Time in seconds')
       if(ibigfigs)
        bigfig;whitefig;boldlines;bigfont;
       end
       [W,fw]=fftrl(pad(w,padpow2(sq)),t);
       plotimage(ones(length(trow),1)*abs(W'),trow,fw);
       title('Wavelet surface');
       xlabel('Frequency in Hz');ylabel('Time in seconds')
       if(ibigfigs)
        bigfig;whitefig;boldlines;bigfont;
       end

       plotimage((ones(length(trow),1)*abs(W')).*alpha.*abs(rqtvs),trow,fw);
       title('Gabor spectrum model');
       xlabel('Frequency in Hz');ylabel('Time in seconds')
       if(ibigfigs)
        bigfig;whitefig;boldlines;bigfont;
       end

    end
end