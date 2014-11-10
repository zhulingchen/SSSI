function [sonicnum,densitynum,waveletnum,densityopt,holealg,a,m,theotype,name]=...
				theodeffini
%
% [sonicnum,densitynum,waveletnum,densityopt,holealg,a,m,theotype,name]=theodeffini
%
% Call this to complete the theogram definition dialog. Return values mean:
%	sonicnum = number of the sonic section selected for theogram
%	densitynum = number of the density section selected for theogram
%	waveletnum = number of the wavelet selected for theogram
%	densityopt = density option selection 1: use Gardners exclusively, 2: use
%		density logs but fill holes with Gardners 3: use density logs exclusively
%	holealg = hole filling algorithm selection 1: constant 2: linear 3: mean
	params=get(gca,'userdata');
	if( params==-1 ) %test for a cancel
		sonicnum=-1;
		return;
	end
	sonicnum=params(1);
	densitynum=params(2);
	waveletnum=params(3);
	densityopt=params(4);
	holealg=params(5);
%		4: layer mean 5: layer trend
%   domain = 1: compute rcs in time  or 2: compute rcs in depth
%	a= scalar coefficient in gardners relation. Will be -1 if densopt=3
%	m= scalar exponent in gardners relation. Will be -1 if densopt=3
%	name = name chosen for the rc section
	a=params(6);
	m=params(7);
	theotype=params(8);
	name=setstr(params(8:length(params)));
