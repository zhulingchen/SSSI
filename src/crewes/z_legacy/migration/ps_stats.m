function ps_stats(nz,j,etm,ets,ttm,tts,rmt,rst)
disp(' ');
disp('  Phase shift migration: v(z)');
disp(['  Depth step ',num2str(j),' of ',nz,]);
%disp(['  Floating point ops: ',num2str(f)]);
disp(['  Iteration time: ',num2str(etm),' min ',num2str(ets),' sec']);
disp(['  Run time so far: ',num2str(ttm),' min ',num2str(tts),' sec']);
disp(['  Approx remaining runtime: ',num2str(rmt),' min ',num2str(round(rst)),' sec']);
