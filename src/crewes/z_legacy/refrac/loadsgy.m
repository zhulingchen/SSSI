[seisp,timep,line_namep,ntrp,nrecp,xsp,ysp,xrp,yrp,offsp,selevsp,relevsp,sdepthsp,cdpsp]=readsegy('soup189.sgy');
[seisgp]=readsegy('sougp189.sgy');
[seisr,timer,line_namer,ntrr,nrecr,xsr,ysr,xrr,yrr,offsr,selevsr,relevsr,sdepthsr,cdpsr]=readsegy('sour189.sgy');
[seisgr]=readsegy('sougr189.sgy');
seisimage(seisgp,timep(1:2001,:))
seisimage(seisgr,timer)
