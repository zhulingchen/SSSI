function [output]=del2_9pt(input,delx);
% DEL2_9PT ... compute the 9 point Laplacian
% 
% [output]=del2_9pt(input,delx)
%
% DEL2_9PT computes the 9 point approximation to the laplacian operator of 
% the two dimensional matrix 'input'.  The horizontal and vertical bin
% spacing of the matrix must be the same.  
%
% input = input matrix
% delx = the horizontal/ vertical bin spacing in consistent units
%
% output = the output matrix
%
% by Carris Youzwishen, April 1999
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

[nrows,ncolumns]=size(input);

input2=zeros(nrows+4,ncolumns+4);
input2(3:nrows+2,3:ncolumns+2)=input;

% create output matrix
output=zeros(nrows,ncolumns);

% for first dimension

  output = (-input2(1:nrows,3:ncolumns+2) + 16*input2(2:nrows+1,3:ncolumns+2) ...
                 -30*input2(3:nrows+2,3:ncolumns+2) + 16*input2(4:nrows+3,3:ncolumns+2) ... 
                 -input2(5:nrows+4,3:ncolumns+2))/(12*delx^2);

% for second dimension 

input2=input2.';
output=output.';

 output = (-input2(1:ncolumns,3:nrows+2) + 16*input2(2:ncolumns+1,3:nrows+2) ...
                 -30*input2(3:ncolumns+2,3:nrows+2) + 16*input2(4:ncolumns+3,3:nrows+2) ... 
                 -input2(5:ncolumns+4,3:nrows+2))/(12*delx^2) + output;

output=output.';
