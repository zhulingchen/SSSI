function loadPMcb
% Call of the ProMax loading function 
[fbtime, fbcoord, shotcoord, nshots, nrecs, uphole, recelev, shotelev] = loadPM;
refdata('clear');
% Save the stuff loaded from the files
refdata('set','fbtime', fbtime);
refdata('set','fbcoord', fbcoord);
refdata('set','shotcoord', shotcoord);
refdata('set','nshots',nshots);
refdata('set','nrecs',nrecs);
refdata('set','recelev',recelev);
refdata('set','shotelev',shotelev);
refdata('set','uphole',uphole);
PMTsetmenus;
