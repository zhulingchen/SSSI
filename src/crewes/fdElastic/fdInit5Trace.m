function [jfr,iweX,isfZ].....
    = fdInit5Trace(iTrAcq,trX,trZ,Dxz,ix1,iz1)
% function [jfr,iweX,isfZ].....
%     = fdInit5Trace(iTrAcq,trX,trZ,Dxz,ix1,iz1)
%Set up for output of surface acquired traces at each time step
%The input parameters are
%iTrAcq  .... Not used (indicator of structured surface)
%trX     .... X-position of the 'well' where trace (time) data will be collected
%trZ     .... Z-depth of the 'line' where trace (time) data will be collected
%Dxz     .... FD sample rate in (metres)
%ix1   ...... X co-ordinate of start of model in arrays
%iz1   ...... Z co-ordinate of start of model in arrays
%The output parameters are
%jfr    ..... The index of the next frame (for movies)
%isfZ    .... From trZ (if no structured surface)
%iweX    .... From trX (if no structured surface)
%
% P.M. Manning, Dec 2011
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

% isurf = 0;
jfr = 0;
% disp(traceFile)
% if ~isempty(traceFile)
%     isurf = round(traceZlevel/Dxz)+iz1;
%     disp(isurf)
%     surfUx = zeros(nxf,nstep);
%     surfUz = zeros(nxf,nstep);
%     surfT = zeros(1,nstep);
%     surfS =  [0 0 shotDepth];
%     surfX = (0:nxf-1)*Dxz;
%     surfZ = ones(1,nxf)*traceZlevel;
% end
if ~isempty(iTrAcq)
    isfZ = round(trZ/Dxz)+iz1;    %Modify if structured surface
end
%     if strcmpi(iTrAcq,'x')
%         nf = nxf;
%         else if strcmpi(iTrAcq,'z')
%                 nf = nzf;
%              end
%     end
if trX > 0
    iweX = round(trX/Dxz)+ix1;
end
% surfUx = zeros(nxf,nstep);
% surfUz = zeros(nxf,nstep);
% wellUx = zeros(nzf,nstep);
% wellUz = zeros(nzf,nstep);
% surfUx = single(zeros(nxf,nstep));
% surfUz = single(zeros(nxf,nstep));
% wellUx = single(zeros(nzf,nstep));
% wellUz = single(zeros(nzf,nstep));
% end