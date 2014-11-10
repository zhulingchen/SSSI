function obj = splitData( obj )
%
% function obj = splitData( obj )
%
% Split data rows based on obj.delimiter using textscan
% We're using %q instead of %s to strip double-quotes from character strings
%

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geoscience of the University of Calgary, Calgary,
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

for i = 1:size(obj.sections,2)
    if obj.isDataSection(i)
        %disp (['data section: ' num2str(i)])
        ca = cellfun(@(X) ...
            textscan(X,'%q','delimiter',obj.delimiter,'MultipleDelimsAsOne', 1),...
            obj.sections{i}{3});
        ca = deblank(ca);
        try
            % re-order to a 2D cell-array with log data in the columns
            ca = [ca{:}]';
        catch ex
            %warning('crewes:las:splitdata',['1 ' ex.message]);
            % displaying exception leads to clutter...
            try
                s = repmat({max(cellfun(@length, ca))},1,length(ca));
                ca = cellfun(@(X, Y) obj.padcell(X, Y), ca, s, 'UniformOutput',false);
                
                ca = [ca{:}]';
            catch ex
               %warning('crewes:las:splitdata',['2 ' ex.message]);
               try
                    ncols = obj.numLogs(i);
                    d=vertcat(ca{:});
                    ca=reshape(d,ncols,length(d)/ncols)';
               catch ex
                    error('crewes:las:splitdata',['3 ' ex.message]);
               end
            end
        end
        obj.sections = {'sectiondata',i,ca};
    end
end %end for
end %end function