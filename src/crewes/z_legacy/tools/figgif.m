% figgif - Save a figure into a GIF file
% figgif by itself saves the current figure into a file 
%        called 'Figure No. 1.gif'
%
% figgif('beep') saves the current figure into a file called 'beep.gif'
%
% figgif('beep',h) saves the figure with the handle 'h' into a file
%        called 'beep.gif'
%
% D. Foltinek, August 1997
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
function figgif(filename, h)
if( nargin < 2 )
   h = gcf;
end
% Bring the figure to the foreground
figure(h);
% Build a string containing the name of the X window
name = get(h, 'name');
if( strcmp(name, '') == 1 )
   figname = sprintf('Figure No. %d ', h );
else
   figname = sprintf('Figure No. %d: %s', h, name );
end
if( nargin < 1 )
   filename = sprintf('%s.gif', figname);
else
   filename = sprintf('%s.gif', filename);
end
% Build the Unix command string to grab the named window, save as GIF
cmdstr = sprintf('xwd -name ''%s'' | xwdtopnm | ppmquant -fs 256 | ppmtogif > ''%s'' ', figname, filename );
unix(cmdstr);
fprintf('Saved figure to %s\n',filename);
