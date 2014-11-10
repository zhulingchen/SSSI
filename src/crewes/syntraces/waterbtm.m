function mtrain=waterbtm(twater,t,rc,flag)
% WATERBTM: compute the sero offset water bottom response
%
% mtrain=waterbtm(twater,t,rc,flag)
% mtrain=waterbtm(twater,t,rc)
% mtrain=waterbtm(twater,t)
%
% WATERBTM computes the impulse response for a zero offset water
%  bottom multiple.
%
% twater= vertical traveltime thickness of the water (2-way)
% t= time coordinate vector. This will determine the sample rate
%    and length of mtrain
% rc= reflection coefficient of the water bottom
% ********** default =.2 *************
% flag= 0 ... produce the combinde response (source and receiver
%             multiples combined)
%       1 ... only the source multiples are desired
% ********** default = 0 ******************
%
% by G.F. Margrave, July 1991
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

 if nargin<4, flag=0; end
 if nargin<3, rc=.2; end
% generate the one way response
 time=twater;
 r=1;
 mtrain=impulse(t,1);
 while time<=max(t)
  tw=near(t,time);
  r=r*(-rc);
  mtrain(tw)=r;
  time=time+twater;
 end
 if flag==0,
   mtrain=conv(mtrain,mtrain);
   mtrain=mtrain(1:length(t));
 end
  
  
 






