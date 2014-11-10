function [delay]=delayt(fbtime,fbcoord,cvpavg,v2rec,shotcoord,nshots,recelev,slim,elim,plust)
delay=NaN*ones(nshots,length(recelev));
for n=1:nshots
    ind=find(plust(1,:)==shotcoord(n));
    ts=plust(2,ind)/2;
    if (length(ts)~=0)
	  validn = find(~isnan(fbtime(n,:)));
	  indl = find(fbcoord(n,validn)<slim & fbcoord(n,validn)<cvpavg(n,1));
	  indr = find(fbcoord(n,validn)>elim & fbcoord(n,validn)>cvpavg(n,2));
	if (length(indl)>0)
	    for h=indl
 	      indcoord(h)=find(recelev(1,:)==fbcoord(n,validn(h)));
 	      delay(n,indcoord(h)) = fbtime(n,validn(h))-ts-(shotcoord(n)-fbcoord(n,validn(h)))/v2rec(1,indcoord(h));
	    end
	end
	if (length(indr)>0)
	    for f=indr
 	      indcoord(f)=find(recelev(1,:)==fbcoord(n,validn(f)));
 	      delay(n,indcoord(f)) = fbtime(n,validn(f))-ts-(fbcoord(n,validn(f))-shotcoord(n))/v2rec(1,indcoord(f));
	    end
	end
    end
end
