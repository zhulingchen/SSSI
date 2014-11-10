% side      = leftside or rightside   1=left  -1=right
% tdiff     = TD curve
% kcoord    = coordinates of the TD curve
% xcoord    = receiver coordinates (all recs in line)
% shotcoord = coordinates of the shot
% tdfold    = input fold array
% tdstack   = array that the TD curves are getting stacked into
% shotnum   = current shot
% tolerance = maximum slope that is considered 'flat' (see fflat)
% maxtolerance = maximum allowable tolerance
% lastvel   = velocity used to flatten the last curve
% spref, epref = start and end indecies of the reference TD curve
% tdref     = reference TD curve
% shotref   = reference shot number
% firstshot = flag ==1 if this is the first shot to be stacked
% secondshot= flag ==1 if this is the second shot to be stacked
function [tdstack, tdfold,tolerance,spref,epref,tdref,shotref,lastvel] = ...
         stackTD(tdiff, kcoord, xcoord, shotcoord, tdfold, tdstack, ... 
                 side, shotnum, tolerance, maxtolerance, lastvel, ...
                 spref, epref, tdref, shotref, firstshot, secondshot)
% Assign a direction flag: coordinates increase with index number or not.
if( kcoord(1) < kcoord(2) )
   rdir = 1;
else
   rdir = -1;
end
if( shotcoord(1) < shotcoord(2) )
   sdir = 1;
else
   sdir = -1;
end
[sp, ep] = fflat(kcoord, tdiff, tolerance, side*rdir);
notfound = (isempty(sp));
if firstshot
   while( notfound & (tolerance < maxtolerance) )
      tolerance = tolerance + 0.002;
      fprintf('stackTD: increasing tolerance to %f\n', tolerance);
      [sp, ep] = fflat(kcoord, tdiff, tolerance, side*rdir);
      notfound = (isempty(sp));
   end
   if( tolerance > maxtolerance )
      disp('Warning: could not find flat portion of TD curve');
      tdstack = NaN * ones(size(tdstack));
      tdfold = NaN * ones(size(tdfold));
      return;
   else
      spref = sp;
      epref = ep;
      tdref = tdiff;
      shotref = shotnum;
   end
else
   % Calculate the distance corresponding to the two TD curves
   sep = shotcoord(shotnum) - (shotcoord(shotref+sdir) - shotcoord(shotref) );
   if notfound
      if secondshot
         while( notfound & tolerance < maxtolerance )
            tolerance = tolerance + 0.002;
            fprintf('stackTD: increasing tolerance to %f\n', tolerance);
            [sp, ep] = fflat(kcoord, tdiff, tolerance, side*rdir);
            notfound = (isempty(sp));
         end
         if( tolerance > maxtolerance )
            disp('Warning: could not find flat portion of TD curve');
            tdstack = NaN * ones(size(tdstack));
            tdfold = NaN * ones(size(tdfold));
            return;
         end
      else
         % Then use the previous velocity to flatten this one.
         fprintf('stackTD: flattening with old velocity: %f\n',lastvel);
         tdiff = tdiff - sep/lastvel;
         plot(kcoord, tdiff,'y');  hold on;
       end
   else
      [tdiff lastvel] = flattenTD(tdiff,tdref,sp,ep,spref,epref,sep);
      fprintf('stackTD: flattening with velocity: %f\n',lastvel);
      plot(kcoord, tdiff,'g');  hold on;
   end
end
if( lastvel ~= NaN & length(lastvel)~=0 )
   for i=1:length(kcoord)
      index = find(xcoord == kcoord(i));
      if(~isempty(index))
         tdstack(index) = tdstack(index) + tdiff(i);
         tdfold(index) = tdfold(index)+1;
      else
         fprintf('stackTD: no match on coordinates to stack shots %d %d\n',...
                 shotnum, i);
      end
   end
end
