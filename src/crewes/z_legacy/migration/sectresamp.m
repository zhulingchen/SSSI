function [seis2,t2] = sectresamp(seis,t,dt2,timesout,flag,fparms)
% SECTRESAMP: runs resamp on each trace in a seismic section
%
%  [seis2,t2] = sectresamp(seis,t,dt2,timesout,flag,fparms)
%
% SECTRESAMP resamples (in time) a seismic section.
% It simply loops over columns and calls resamp for each trace.
%
% seis ... input section of size nsamp x ntr. That is one trace per
%	column.
% t ... nsamp long time coordinate vector for seis
% dt2= desired output sample size (in seconds). May be any floating point 
%           number. 
%   timesout = start and end time desired for the output trace: [tmin tmax]
%  	If dt2 does not divide evenly into this interval then the first
%  	These times need not fall on the input time grid represented by t
%  	sample will be at tmin while the last will be just less than tmax
%         ********* default [t(1) t(length(t))] **************
%   flag= antialias filter option needed if dt2> t(2)-t(1) 
%         0 -> a zero phase antialias filter is used
%         1 -> a minimum phase antialias filter is used
%         ***************** default = 1 *****************
%   fparms= 2 element vector of antialias filter parameters
%         fparms(1)= 3db down point of high frequency rolloff
%         expressed as a fraction of the new nyquist
%         fparms(2)= attenuation desired at the new Nyquist in db
%   *************** default = [.6,120] *******************
%
% seis2 ... output section of size nsamp x ntr.
% t2 ... output time coordinate vector
%
% G.F. Margrave, CREWES Project, University of Calgary, 1996
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


[nsamp,ntr]=size(seis);

if(nargin < 4)
	for k=1:ntr
		[tmp,t2] = resamp( seis(:,k),t,dt2);
		if(k==1)
			seis2=zeros(length(t2),ntr);
		end
		seis2(:,k)=tmp;
		if(rem(k,50)==0)
			disp(['finished ' num2str(k) ' traces'])
		end
	end
elseif( nargin < 5)
	for k=1:ntr
		[tmp,t2]  = resamp( seis(:,k),t,dt2,timesout);
		if(k==1)
			seis2=zeros(length(t2),ntr);
		end
		seis2(:,k)=tmp;
		if(rem(k,50)==0)
			disp(['finished ' num2str(k) ' traces'])
		end
	end
elseif( nargin < 6)
	for k=1:ntr
		[tmp,t2]  = resamp( seis(:,k),t,dt2,timesout,flag);
		if(k==1)
			seis2=zeros(length(t2),ntr);
		end
		seis2(:,k)=tmp;
		if(rem(k,50)==0)
			disp(['finished ' num2str(k) ' traces'])
		end
	end
else
	for k=1:ntr
		[tmp,t2]  = resamp( seis(:,k),t,dt2,timesout,flag,fparms);
		if(k==1)
			seis2=zeros(length(t2),ntr);
		end
		seis2(:,k)=tmp;
		if(rem(k,50)==0)
			disp(['finished ' num2str(k) ' traces'])
		end
	end
end
