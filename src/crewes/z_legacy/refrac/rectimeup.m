% nshots = number of shots 
% nrecs = number of receivers 
% reclocation = receiver location matrix (nshots,nrecs)
% shotpick = fb pick matrix (nshots, nrecs)
% r = receiver location vector (nrecs)
%
% Output: diffmat = matrix of reciprocal shot fb pick differences
function diffmat = rectimeup(nshots, nrecs, reclocation, pickuphole, r)
diffmat(nshots-1, nshots) = NaN;
for i=1:nshots-1;       % i is the first shot (t1)
  for j=i+1:nshots;    % j is the 2nd shot (t2)
    t1 = NaN;
    t2 = NaN;
    if( (r(j)>min(reclocation(i,:))) & (r(j)<max(reclocation(i,:))) )
       t1=interp1(reclocation(i,:),pickuphole(i,:),r(j),'linear');
    end
    if( r(i) > min(reclocation(j,:)) & r(i) < max(reclocation(j,:)) )
       t2=interp1(reclocation(j,:),pickuphole(j,:),r(i),'linear');
    end
    diffmat(i,j) = t2-t1;
  end
end
