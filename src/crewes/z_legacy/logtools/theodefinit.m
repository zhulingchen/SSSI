function theodefinit(transfer,hmasterfig,sonics,densities,wavelets,...
                     densopts,holealg,a,m,msg,sonicnum,densitynum,wletnum,...
                     sectionopt,typeopt,nameopt,name)
%
% theodefinit(transfer,hmasterfig,sonics,densities,wavelets,...
%             densopts,holealg,a,m,msg,sonicnum,densitynum,wletnum,...
%             sectionopt,typeopt,nameopt,name)
%
% THEODEFINIT is called by LOGEDIT to initiate a dialog  to define the
% computation of reflection coefficient sections from sonic log(s)
% and, optionally, density log(s).
%
% transfer ... transfer command to be called when the user terminates
%              the dialog
% hmasterfig ... handle of the masterfigure (usually a LOGEDIT window) in
%                control of the dialog
% sonics ... list of names of possible sonic logs to be used to compute the
%            theogram. Should be a string row vector with individual 
%            names separated by '|' such as: 'fred|sam|billy|wilma'
% densities ... list of names of possible density logs to be used to compute
%               the theogram. Should be a string row vector with individual 
%               names separated by '|' such as: 'fred|sam|billy|wilma'
% wavelets ... list of names of possible wavelets to be used to create the
%              theogram. Should be a string row vector with individual names
%              separated by '|' such as: 'fred|sam|billy|wilma'
%              If null ('') then no wavlet information will be requested.
% densopts ... flag giving the default for the density option:
%              0= use constant density
%              1= use Gardners relation exclusivly; 
%              2= use density section where defined, fill in holes with 
%                 Garnders relation after hol filling on sonic section 
%              3= Use density section exclusively. Fill holes on it
%		  independently of the sonic section.
% holealg ... flag signifying the preferred default for the holefilling 
%             algorithm: 
%             1=constant; 2=linear; 3=mean; 4=layermean; 5=layertrend;
%             See FILLHOLES for more information
% Gardners relation assumes density= a*(vins).^m
% a ... default scalar multiplier in Gardner's relation
% m ... default scalar exponent in Gardner's relation
% msg ... a message to be displayed at the top of the dialog
% **************** default '' ****************
% sonicnum ... integer denoting the initial selection in sonics for the 
%              sonic log
% *************** default 1 ***************
% densitynum ... integer denoting the initial selection in densities for 
%                the density log
% *************** default 1 ***************
% wletnum ... integer denoting the initial selection in wavelets for the 
%             wavelet
% *************** default 1 ***************
% sectionopt ... if 0, then log labels refer to single logs, if 1 they 
%                refer to sections
% *************** default 0 ***************
% typeopt ... if 0, then no theogram type choice is given
%             if 1, then a theogram with primaries only is requested.
%             if 2, then theoretical 1-d seismogram with attenuated
%             primaries plus multiples is requested.
% *************** default 0 ***************
% nameopt ... if 0, then no name is requested for the theogram.
%		if 1, then a name must be provided
% *************** default 0 ***************
% name ... string giving the initial name choice
% *************** default is automatically generated name *********
%
% G.F. Margrave, November 1994
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
% theodefinit(transfer,hmasterfig,sonics,densities,wavelets,...
%			densopts,holealg,a,m,msg,sonicnum,densitynum,wletnum,...
%			sectionopt,typeopt,nameopt,name)
if(nargin<17)
	name=0;
end
if(nargin<16)
	nameopt=0;
end
if(nargin<15)
	typeopt=0;
end
if(nargin<14)
	sectionopt=0;
end
if(nargin<13)
	wletnum=1;
	if( isempty(strcmp(wavelets)) )
		wletnum=0;
	end
	if(isempty(wavelets))
		wletnum=0;
	end
end
if(nargin<12)
	densitynum=2;
end
if(nargin<11)
	sonicnum=1;
end
if(nargin<10)
	msg='';
end
% pack the information into the current axes userdata
hax=get(hmasterfig,'currentaxes');
figure(hmasterfig);
set(hax,'userdata',[abs(transfer) nan abs(msg) nan abs(sonics) nan...
	abs(densities) nan abs(wavelets) nan hmasterfig densopts holealg...
	a m sonicnum densitynum wletnum sectionopt typeopt nameopt abs(name)]);
theodef('init')
