function selboxinit(action,button)
% SELBOXINIT: initialize selection box drawing
%
%  selboxinit(transfer,button)
%  selboxinit(transfer)
%  selboxinit
%
% initialize the current figure and axes for a selectionbox. See help for
% SELBOX.
% transfer ... MATLAB string containing an executable command which is 
%		called at the completion of the selection box
%	********** default = '' ************
% button ... mouse button number to respond to the selection box. If the
%		user presses a different button than this, no action is taken
%	********** default = 1 ************
%
% G.F. Margrave, CREWES, November 1995
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

%the following function pragma is needed for the compiler
%#function sca

if(nargin<2)
		button=1;
	end
	if( button== 1)
		set(gcf,'windowbuttondownfcn','sca;selbox(''init'',1)');
		if( nargin < 1)
			set(gcf,'windowbuttonupfcn','selbox(''fini'',1)');
		else
			set(gcf,'windowbuttonupfcn',...
				['selbox(''fini'',1);' action]);
		end
	elseif(button==2)
		set(gcf,'windowbuttondownfcn','sca;selbox(''init'',2)');
		if( nargin < 1)
			set(gcf,'windowbuttonupfcn','selbox(''fini'',2)');
		else
			set(gcf,'windowbuttonupfcn',...
				['selbox(''fini'',2);' action]);
		end
	elseif(button==3)
		set(gcf,'windowbuttondownfcn','sca;selbox(''init'',3)');
		if( nargin < 1)
			set(gcf,'windowbuttonupfcn','selbox(''fini'',3)');
		else
			set(gcf,'windowbuttonupfcn',...
				['selbox(''fini'',3);' action]);
		end
	end
