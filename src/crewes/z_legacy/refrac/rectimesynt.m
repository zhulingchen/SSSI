% nshots = number of shots 
% nrecs = number of receivers 
% fbcoord = receiver location matrix (nshots,nrecs)
% fbtime = fb pick matrix (nshots, nrecs)
% shotcoord = shot location vector (nshots)
%
% Output: diffmat = matrix of reciprocal shot fb pick differences
function diffmat = rectimesynt(nshots, nrecs, fbcoord, fbtime, shotcoord)
diffmat(nshots-1, nshots) = NaN;
for i=1:nshots-1;       % i is the first shot (t1)
  for j=i+1:nshots;    % j is the 2nd shot (t2)
    t1 = NaN;
    t2 = NaN;
    if( (shotcoord(j)>min(fbcoord(i,:))) & (shotcoord(j)<max(fbcoord(i,:))) )
       t1=interp1(fbcoord(i,:),fbtime(i,:),shotcoord(j),'linear');
    end
    if( shotcoord(i) > min(fbcoord(j,:)) & shotcoord(i) < max(fbcoord(j,:)) )
       t2=interp1(fbcoord(j,:),fbtime(j,:),shotcoord(i),'linear');
    end
    diffmat(i,j) = t2-t1;
  end
end
