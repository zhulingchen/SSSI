This directory contains the matlab scripts that can
read and write data into SAC format.
Among them, fget_sac.m, sachdr.m is written by me,
sac.m is written by Dr. Xianglei Huang at umich.
I forgot where I got the rest from (maybe Jeff McGuire at WHOI?!).

Program		Description
fget_sac.m	main code to load sac data into matlab
sachdr.m	subroutine to convert the sachdr into readble struture array
sac.m		subroutine to load binary sac data

newSacHeader.m	subroutine to generate a sac hdr
rdSacHead.m	read header of SAC format data
rdSac.m		read SAC format data
wtSac.m		write SAC format data

sacfft.m	a subroutine to compute FFT for SAC data

Example:

1. load SAC data N.MYJH.Z.sac using fget_sac.m
Command: [Ztime,Zdata,ZSAChdr] = fget_sac('N.MYJH.Z.sac');
Note: type ZSAChdr to take a look at the headers, to use info. p 
      arrival in the header, type ZSAChdr.times.a.

2. write ascii data (MYJH.dat) into SAC in matlab using wtSac.m
Command: 
load MYJH.dat; % load the ascii data
N = length(MYJH); % get the total length
dt = 0.01; % sampling rate
tstart = 0; % starting time
MYJH_hd = newSacHeader(N,dt,tstart);
MYJH_sacfile = 'MYJH.sac';
wtSac(MYJH_sacfile,MYJH_hd,MYJH);

Note: 
1. If you want to add info. into headers, you have to do it
      by hand, check the binary header format at
      http://www.iris.edu/manuals/sac/SAC_Manuals/FileFormatPt1.html 

2. Other sac matlab program can be found at
http://www.aeic.alaska.edu/input/mthorne/software/index.html

Last updated by zpeng, Sun Nov  5 17:52:50 EST 2006
