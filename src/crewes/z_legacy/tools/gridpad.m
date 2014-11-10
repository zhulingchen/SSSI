function grout = gridpad( grin, nrows, ncols, padval)
% grout = gridpad( grin, nrows, ncols, padval)
% grout = gridpad( grin, nrows, ncols)
%  
% pad (or truncate) a grid to a new size
%
% grin = input grid
% nrow = number of rows on output
% ncols = number of columns on output
% padval = value to pad with
%		*********** default 0.0 ***********
% grout = output grid
%
% by G.F. Margrave October 1993
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
if( nargin < 4) padval = 0.0; end
[rin colin]=size(grin);
% pad the rows
	if( nrows <= rin ) % truncate
		grout = grin(1:nrows,:);
	else %pad
		grout = [grin; padval*ones(nrows-rin,colin)];
	end
	
% pad the columns
	if( ncols <= colin ) %truncate
		grout = grout(:,1:ncols);
	else % pad
		grout = [grout padval*ones(nrows,ncols-colin)];
	end
