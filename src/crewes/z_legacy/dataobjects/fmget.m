function val=fmget(fmat,action)
%
% val=fmget(fmat,action)
%
% FMGET is used to retrieve information from a flexi-matrix (see FMSET for a 
% definition. There are three possible values for action which result in three
% different return values:
%
% action = 'x' ... returns the x coordinates (positions of the columns) 
%                  of the matrix
% action = 'y' ... returns the y coordinates (positions of the rows) of 
%                  the matrix
% action = 'mat' ... returns the matrix
%
% Note that, due to the definition of a flexi-matrix, y values will be a 
% regular ordered set while the x values may be completely irregular 
% and in any order.  Also, the matrix will most likely contain a number of
% nan's in each column as these are used to pad the columns to equal length
%
% by G.F. Margrave
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

[nr,nc]=size(fmat);

if(nr*nc==0)
	val=[];
	return;
end

if(strcmp(action,'mat'))
	val=fmat(2:nr,2:nc);
	return;
end

if(strcmp(action,'x'))
	val=fmat(1,2:nc);
	return;
end

if(strcmp(action,'y'))
	val=fmat(2:nr,1);
	return;
end
	
