function newlineinit(transfer)
% newlineinit(transfer)
%
% NEWLINEINIT initiates the creating of new lines in a figure window using
% NEWLINE.  The sole argument, transfer, should be a MATLAB string with a
% syntactically correct MATLAB command to be executed when a new line has
% been created. For example if the function my_draw is calling newlineinit
% to initiate drawing a line, then when NEWLINE finishes a line it would be
% desirable to have it call my_draw with an appropriate argument such as the
% string 'linedone'. This could be accomplished by:
% newlineinit('my_draw(''line_done'')')
%
% The finished line and its graphic handle will be found in the userdata
% of the current axes (See help for newlinefini for details). 
% Therefore it is important that my_draw ensure there is
% no important information stored there when newlineinit. To retrieve the
% newline coordinates, my_draw must call newlinefini because some additional
% cleanup must be done.  NEWLINE will terminate by calling my_draw everytime
% it is signalled that a new line with a unique graphics handle has been
% created. This occurs whenever:
%	1) A new line has been drawn with mouse button 1 and terminated by
%	   pressing button 2 or 3
%	2) A new line has been created by duplicating an existing line
%	3) A line has been deleted. In this case its handle is returned as
%	   negative
%	4) A deleted line has been restored with newline('undo')
%
% While NEWLINE has its own undo, it is often more convenient to do it
% yourself.
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
	if(nargin<1)
		transfer=[];
	end
	newline('init',transfer);
