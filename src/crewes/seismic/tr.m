function trout=tr(trin,flag)
% TR Reverses the a time series
% trout=tr(trin,flag)
%
% TR returns the time reverse of trin.
% This is accomplished such that the sample at length(trin)/2
% remains fixed. For length(trin) even, the last sample in trin
% is discarded
%
% trin= input trace
% flag= 0 ... keep the center sample (at length(trin)/2) fixed.
%             For length(trin) even, the last sample is discarded.
%     = 1 ... don't worry about the center sample, simply reverse
%             the time series.
%     ************ default= 0 *************
% trout= output trace
%
% See also FLIPDIM
%
% by G.F. Margrave, June 1991
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

if nargin<2
    flag=0; 
end

%trout=trin(length(trin):-1:1);

% if flag==0,
%     if rem(length(trin),2)==0,
%         trout=[trout(2:length(trin)) 0];
%     end
% end

trout=ifft(conj(fft(trin)));
