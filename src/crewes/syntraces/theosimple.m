function [theo,tlog,rcs,pm,p]=theosimple(sonic,dens,z,fmult,wlet,tw)
% [theo,tlog,rcs,pm,p]=theosimple(sonic,dens,z,fmult,wlet,tw)
%
% THEOSIMPLE computes a 1-D seismic model (Theogram) from well
% log information. It uses the sonic log to build a time-depth
% curve and then calls THEOGRAM. (See THEOGRAM help for more info)
%
% sonic ... vector of sonic log samples. Must be in typical sonic
%	units (e.g. micro-sec/meter).
% dens ... vector of same length as sonic giving density values.
%	For a constant density result, provide this as ones(size(sonic)).
% z ... vector of same size as sonic giving the depths of the log.
% fmult ... multiple flag
%            0 -> produce a primaries only convolutional theogram
%            1 -> produce a primaries (with attenuation losses) plus
%                  multiples theogram
% wlet ... vector of wavelet samples
%   Wavelets can be generated using any of WAVEMIN, WAVEZ, WAVEVIB,
%   RICKER, and FILTF in the Seismic_Toolbox
% tw ...  vector of time coordinates for the wavelet. Note that
%         a non-causal wavelet should have tw(1)<0. The convolution
%         will be such that the wavelet sample at tw==0 will be
%         superimposed on the reflectivity spike
%
% theo ... n length vector containing the theogram samples.
%         (n is determined by the algorithm).
% tlog ... nlength vector containing the times for the theogram
%	rcs ... n length vector of computed reflection coefficients in time
%  pm ... nlength vector of attenuated primaries plus multiples. 
%         If fmult == 0, this is identical to rcs.
%	p ... nlength vector with attenuated primary rcs. If fmult==0, this
%         is identical to rcs. If fmult==1, then the multiple content can
%         be obtained as: pm-p
%
% G.F. Margrave, October 1994
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
% a bit of error checking
	if(length(sonic)~=length(dens))
		error('Sonic and density logs must be the same size');
	end
	if( length(sonic)~=length(z) )
		error('Sonic and Z vector must be the same size');
	end
	if( length(wlet)~=length(tw) )
		error('Wavelet and TW must be the same size');
	end
%obtain a time-depth function
nlegs=length(sonic)/100;
if(nlegs>150) nlegs=150; end
if( nlegs< 10) nlegs=10; end
[tz,zt]=sonic2tz(sonic,z,nlegs);
%convert sonic to vins
sonic=1.e06 ./ sonic;
%call theogram
[theo,tlog,rcs,pm,p]=theogram(sonic,dens,z,wlet,tw,[zt tz],0,fmult);
