function [sonicnum,densitynum,densityopt,holealg,theopt,a,m,name]=...
				theodefinefini
%
% [sonicnum,densitynum,densityopt,holealg,theopt,a,m,name]=theodefinefini
%
% Call this to complete the theogram definition dialog. Return values mean:
%	sonicnum = number of the section selected for sonic logs
%	densitynum = number of the section selected for density logs
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
	densityopt=params(3);
	holealg=params(4);
	theopt=params(5);
%		4: layer mean 5: layer trend
%   theopt = 1: primaries only  or 2: primaries plus multiples
%	a= scalar coefficient in gardners relation. Will be -1 if densopt=3
%	m= scalar exponent in gardners relation. Will be -1 if densopt=3
%	name = name chosen for the rc section
	a=params(6);
	m=params(7);
	name=setstr(params(8:length(params)));
