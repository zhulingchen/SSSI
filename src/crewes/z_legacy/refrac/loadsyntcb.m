function loadsyntcb
% Call the synthetic loading function
[fbtime, fbcoord, shotcoord, nshots, nrecs, recelev, shotelev, uphole] = loadsynt;
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
% Update menus
PMTsetmenus;
