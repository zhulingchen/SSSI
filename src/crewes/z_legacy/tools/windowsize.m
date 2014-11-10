function windowsize(w,h,hfig)
%
% windowsize(w,h,hfig)
% windowsize(w,h)
% windowsize
%
% WINDOWSIZE sets the size of a figure window to the specified 
% w and h, in pixels. Windowsize by itself gives the current 
% window size in pixels.
% 
% w	= width in pixels
% h	= height in pixels
% hfig	= figure handle
% ================== Default = current figure =================
%
% Darren Foltinek, CREWES Project, 1996
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
if( nargin < 3 )
   hfig = gcf;
end
pos = get( hfig, 'position' );
oldx = pos(1);
oldy = pos(2);
if( nargin < 2 )
   disp( sprintf('Window is currently %d x %d pixels\n', pos(3), pos(4) ) );
   disp('Usage: windowsize( X_size, Y_size)');
   return;
end
set( hfig, 'position', [oldx oldy w h]);
