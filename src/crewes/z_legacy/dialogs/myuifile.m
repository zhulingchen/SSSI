function [file,path]=myuifile(hmasterfig,filtspec,dtitle,uitype,varargin)
% [file,path]=myuifile(hmasterfig,filtspec,dtitle,type)
% [file,path]=myuifile(hmasterfig,filtspec,dtitle,type,defaultfile)
%
% Provides a sane interface to uigetfile and uiputfile in that the
% position of the dialog box is adjusted such that its upper left
% corner is aligned with the upper left corner of the master figure.
% 
% hmasterfig     handle of the master figure
% filtspec       string specifying the filter specification in the file
%		         dialog. Will be the first argument to uigetfile or uiputfile.
% dtitle         string specifying the dialog title
% type           must be either 'get' or 'put' and determines whether uigetfile
%		         or uiputfile is called
% defaultfile    [optional] if present, the default filename for the 
%                get/put diaglog.  May also contain a directory name.
%
% A cancel is indicated by a value of 0 for 'file'
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

%determine the position of the calling figure
u=get(hmasterfig,'units');
set(hmasterfig,'units','pixels');
% pos=get(hmasterfig,'position');
set(hmasterfig,'units',u);

u=get(0,'units');
set(0,'units','pixels');
% scsz=get(0,'screensize');
set(0,'units',u);

% Note: due to a compiler bug with uigetfile, we ignore these x,y
% when a default file is provided.
% x=pos(1);
% y=scsz(4)-pos(2)-pos(4);

%assume failure
file=0;
path=0;

%branch on type
if(strcmp(uitype,'get'))
    if (nargin >= 5 && ~isempty(varargin{1}))
        defaultfile = varargin{1};
        [file,path]=uigetfile(filtspec,dtitle,defaultfile);
    else
        [file,path]=uigetfile(filtspec,dtitle);
    end
elseif(strcmp(uitype,'put'))
    if (nargin >= 5 && ~isempty(varargin{1}))
        defaultfile = varargin{1};
        [file,path]=uiputfile(filtspec,dtitle,defaultfile);
    else
        [file,path]=uiputfile(filtspec,dtitle);
    end
end
    
end %end function
