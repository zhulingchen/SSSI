function h=drawpick(x0,t0,dtdx,dx,figno,col,lw)
% DRAWPICK: draw a line on top of an image (e.g. seismic)
%
% drawpick(x0,t0,dtdx,figno,dx,col,lw)
%
% This function makes it easy to draw a "pick" on top of a seismic
% image. The pick is defined by three numbers (x0,t0,dtdx) that give
% the center of the pick and the timedip. The pick is drawn with the
% specified timedip centered at (x0,t0).
%
% (x0,t0) ... position of the center of the line
% dtdx ... timedip of the line
%	******** default 0 ********
% dx ... horizontal dimension of the line
%	******** default width-of-axes/20 *******
% figno ... figure number to draw in
%       ******** default gcf ********
% col ... color to draw with
%	******** default 'r' *******
% lw ... line thickness to draw with
%	******** lw = 2 **********
%
% h ... handle of the drawn line
%
% G.F. Margrave, CREWES, July 2000
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



if(nargin<3) dtdx=0; end
if(nargin<4)
	hax=get(figno,'currentaxes');
	xlims=get(hax,'xlims');
	xw=diff(xlims)/20;
end
if(nargin<5) figno=gcf; end
if(nargin<6) col='r'; end
if(nargin<7) lw=2; end
	
x1=x0-dx/2;
x2=x0+dx/2;
t1=t0+dtdx*(x1-x0);
t2=t0+dtdx*(x2-x0);
figure(figno);

h=line([x1 x2],[t1 t2],[1 1],'color',col,'linewidth',lw)
