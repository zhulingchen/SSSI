function fdataout=ps_rezero(fdata,f,dx,tmax)
% PS_REZERO ... rezero the temporal zero pad in phase shift extrapolation
%
% fdataout=ps_rezero(fdata,f,dx,tmax)
%
% Designed for use with ips and pspi_mig etc.
%
% fdata ... fk spectrum of data with k axis wrapped and positive f's only
% f ... f coordinate vector for fdata
% dx ... spatial sample size for fdata
% tmax ... time (in seconds) beyond which the data will be set to zero
% fdataout ... fkspectrum of data with the zero pad re-zero'd
% 

%ok, we must do an ifft over f, re-zero, and then an fft over t

%figure out tmax and dt
[r,c]=size(fdata);
df=f(2)-f(1);
Tmax=1/df; %true maximum time including any zero pad
fmax=f(end);
dt=.008;
fnyq=.5/dt;
while(fnyq<fmax)
    dt=dt/2;
    fnyq=.5/dt;
end
Tmax=Tmax-dt;%necessary fiddle
n=round(tmax/dt)+1;
%make f and k axes
fnew=0:df:fnyq;
kx=fftshift(1/2/c/dx*[-c:2:c-2]);
data=zeros(length(fnew),c);
indf=near(fnew,f(1),f(end));
data(indf,:)=fdata;
data(end,:)=0;
%inverse fk transform
[tdata,t,x]=ifktran(data,fnew,kx,0,0);
indt=near(t,tmax,Tmax);
tdata(indt,:)=0;
%forward fk transform
data=fktran(tdata,t,x,0,0,0,0);
fdataout=data(indf,:);



