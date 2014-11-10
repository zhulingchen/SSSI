% flattenTD - flattens the TD curve (td) to the reference TD curve
% (tdref) and returns the flattened curve and velocity corresponding
% to the timeshift required to flatten them.
%
% Input:
%   td     - TD curve to flatten
%   tdref  - reference curve.  TD is flattened to tdref
%   sp, ep - start and end index of the 'flat' portion of TD
%   spref, epref - start and end index of the flat portion of TDref
%   x      - seperation between the shots for TD and TDref
function [newtd, vel] = flattenTD(td, tdref, sp, ep, spref, epref, x)
if( sp<ep )
   dir = 1;
   p1 = max(sp, spref);
   p2 = min(ep, epref);
   ok = (p1<p2);
else
   dir = -1;
   p1 = min(sp, spref);
   p2 = max(ep, epref);
   ok = (p1>p2);
end
% The 'ok' flag is true if the two curves overlap
if( ok )
   fprintf(1, 'flatten: p1:%d  p2:%d\n', p1, p2);
   timeshift = mean(td(p1:dir:p2) - tdref(p1:dir:p2));
   vel = x / timeshift;
   newtd = td - timeshift;
   fprintf(1,'flatten: timeshift %f  vel: %f\n', timeshift, vel);
else
   fprintf(1,'flattenTD: curves did not overlap.\n');
end
