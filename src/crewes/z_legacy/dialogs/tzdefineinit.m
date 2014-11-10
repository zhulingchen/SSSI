function tzdefineinit(transfer,hmasterfig,sections,holealg,numlegs,...
			vo,ax,az,isec,secmeans)
% tzdefineinit(transfer,hmasterfig,sections,holealg,numlegs,vo,ax,az,isec)
%
% TZDEFINEINIT is called by LOGSEC to initiate a dialog  to define the
% computation of T-Z functions from sonic logs.
%
% transfer ... transfer command to be called when the user terminates
%              the dialog
% hmasterfig ... handle of the masterfigure (usually a LOGSEC window) in
%                control of the dialog
% sections ... list of names of possible sonic log sections to be integrated
%              to compute the T-Z functions. Should be a string row vector 
%              with individual names separated by '|' as is returned by 
%              objget(anyobject, 'fieldnames')
% holealg ... flag signifying the preferred default for the holefilling 
%         algorithm: 1=constant; 2=linear; 3=mean; 4=layermean; 5=layertrend;
%         See FILLHOLES for more information
% ********* default = 2 *********
%
% numlegs ... the default number of legs for the t-z function approximation
% ********* default = 30 *********
%
% Velocity above the upper LOGSEC layer is assumed specified in the form
% v(x,z) =vo +alphax*x + alphaz*z
% where alphax is an accelerator for x and alphaz is an accelerator for z.
% vo ... constant upper velocity. A value of inf will cause the upper 
%	LOGSEC horizon to be at time 0
%   ******** default inf **********
% alphax ... x accelarator
%   ********* default 0.0 **********
%
% alphaz ... z accelerator 
%   ********** default 0.0 **********
% isec ... number of the default choice for the section name
%   ********* default =1 *******
%
% G.F. Margrave, March 1994
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
if( nargin < 10), secmeans=[.6*ones(length(sections)),zeros(length(sections))]; end
if( nargin < 9), isec=1; end
if(nargin < 8)
	az=0.;
end
if( nargin < 7), ax=0.; end
if( nargin < 6), vo=inf; end
if( nargin < 5), numlegs=30; end
if(nargin < 4), holealg=2; end
% pack the information into the current axes userdata
hax=get(hmasterfig,'currentaxes');
set(hax,'userdata',[abs(transfer) nan abs(sections) nan ...
	hmasterfig holealg numlegs vo ax az isec secmeans]);
tzdefine('init')
