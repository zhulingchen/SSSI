% Function fflat - Finds the portion of a curve, starting at one
% end, which is flat to within a slope of 'tolerance'
%
% Usage:   [p1 p2] = fflat(x,y,tolerance,beginning);
%
%      x, y - data defining the curve
% tolerance - the maximum slope of the flat portion
% beginning - 1==start at minimum coordinate  -1==start at maximum coord
% 
% p1, p2 - indexes of the flat portion of the curve
% If a straight segment cannot be found, then p1 and p2 
% are both returned empty ([])
function [p1, p2] = fflat(x, y, tolerance, beginning)
curvestart = 1;
curveend = length(y);
if( beginning == 1 )
   p1 = curvestart;
   p2 = curveend;
   dir = 1;
else
   p1 = curveend;
   p2 = curvestart;
   dir = -1;
end
% If the minimum x-coordinate is at the start of the x array (x(1))
% then the coordinates increase with the index numbers
xminindex = find(x == min(x));
if( xminindex == 1)
   reverse = 0;
else
   reverse = 1;
end
% If the coordinate direction is reversed, swap p1 and p2,
% and set the search direction to be the other way
% if( reverse )
%    tmp = p1;
%   p1 = p2;
%   p2 = tmp;
%   dir = -dir;
% end
p2max = curveend-1;
p2min = curvestart+1;
dist = dir*length(y);
done = 0;
while ~done
   lastp2 = p2;
   [p tmp] = polyfit(x(p1:dir:p2), y(p1:dir:p2), 1);
   slope = p(1);
   dist = dir*floor(abs(dist)/2);
%   fprintf(1, 'slope: %f  p2: %d dist: %d\n', slope, p2, dist);
   if( abs(slope) < tolerance )
      p2 = floor(p2 + dist);
   else
      p2 = floor(p2 - dist);
   end 
   p2 = min(p2, p2max);
   p2 = max(p2, 2);
   done = (abs(dist)==1 | p2==p2max | p2==p2min );
end
% Now we have a close answer.  
% Do a linear search until we hit the tolerance limit.
if abs(slope) < tolerance 
   sdir = dir;
%   fprintf('Doing linear search to increase length\n');
else
   sdir = -dir;
%   fprintf('Doing linear search to decrease error\n');
end
done = 0;
while ~done
   lastp2 = p2;
   lastslope = slope;
   p2 = min(p2+sdir, p2max);
   p2 = max(p2,p2min);
   [p tmp] = polyfit(x(p1:dir:p2), y(p1:dir:p2), 1);
   slope = p(1);
%   fprintf(1, 'slope: %f  p2: %d dist: %d\n', slope, p2, dist);
   if(sdir == dir)
      done = (abs(slope)>tolerance);
   elseif(sdir == -dir)
      done = (abs(slope)<tolerance);
   end
   done = done | (p2==p2max) | (p2==p2min);
end
% The linear search will go one point too far.  Backup if needed.
if( abs(slope) > tolerance )
   p2 = lastp2;
   slope = lastslope;
end
fprintf(1, 'final slope: %f  p2: %d dist: %d\n', slope, p2, dist);
% Define the 'flat spot not found' condition.
if( (abs(slope) > tolerance) | (abs(p1-p2) < 5) )
%   figure; 
%   plot(x,y, 'g');
   p1 = [];
   p2 = [];
else
%   line([x(p1) x(p2)], [y(p1) y(p2)]);
end
