function editlinesinit(anchorinfo)
% EDITLINESINIT: called to initiate editing with EDITLINES.
%
% To use EDITLINES in a tool of your own crafting, it is usually necessary 
% to have the graphics handles of the lines to be edited if precise control
% over the editing process is to be maintained. Suppose that there are three
% (x,y) curves on display in the current axes of the current figure and let
% their handles be hline1, hline2, and hline3 and let their coordinates be
% given by vectors: x1,y1,x2,y2,x3,y3.  Then the following sequence
% initiates  editing, terminates editing, and retrieves the edited (x,y)'s:
%
% editlinesinit; .... Begin editing. All three curves are available for
%                     editing. Editing occurs for an indefinite period of time
% editlinesfini; .... Terminate the editing session
% x1=get(hline1,'xdata'); ... get line 1's x coodinates after editing
% y1=get(hline1,'ydata'); ... get line 1's y coodinates after editing
%
% By supplying an argument to editlinesinit, you can control which curves
% may be edited and if there are any initial anchor points on the curves.
% (See anchor points below). (Note that editlinesinit with no argument means
% everything is editable).  For example:
%
% editlinesinit([hline1 0 hline3 0]) ... line 2 may not be edited. Lines
%                                        1&3 have no anchor points
% editlinesinit([hline1 1 x1(10) y1(10) hline2 2 x2(5) y2(5) x2(15) ...
%               y2(15) hline3 0])
%	... this means all three lines may be edited. Line one has a single
%	    anchor at the point x1(10) y1(10); line2 has 2 anchors at the
%	    fifth and 15th points; while line 3 has no anchors.
%
% The argument to editlinesinit is called the anchor vector and gives
% quite general control over editing. editlinesfini returns an updated
% anchor vector in the same format as the input but with any changes which
% the user has made. Thus it can be saved and supplied to editlinesinit to
% resume editing in the same state. 
%
% There are also several different modes that editlines can be run in.
% Modes are activated at anytime in a editing session by an apporpriate 
% call to editlines:
%
% editlines('xonly') ... changes cannot be made to the y coordinates of 
%                        any point.
% editlines('yonly') ... changes cannot be made to the x coordinates of 
%                        any point.
% editlines('xandy') ... frees both x and y for changes
% editlines('linkmode') ... toggles between link mode and normal mode
% editlines('dragmode') ... toggles between constant drag and elastic drag
% editlines('constantdrag') ... sets drag mode to constant
% editlines('elasticdrag') ... sets drag mode to elastic
%
% Note that these calls are on-off toggles. That is if editing is already 
% in 'xonly' mode then: editlines('xonly') will free up the x coordinates 
% for editing.   Also, if both 'xonly' and 'yonly' are on then no changes 
% can be made.
%
% Lastly, editlines('undo') will undo changes made since the last mouse down.
%
% G.F. Margrave December 1993
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

	set(gcf,'windowbuttondownfcn','editlines(''buttondown'')');
	set(gcf,'windowbuttonupfcn','');
	set(gcf,'windowbuttonmotionfcn','');
	
	if(nargin>0)
		set(gca,'userdata',anchorinfo);
	else
		set(gca,'userdata',[]);
	end
	
	editlines('init');


