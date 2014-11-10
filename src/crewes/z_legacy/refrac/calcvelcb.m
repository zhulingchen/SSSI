function calcvelcb
% Function calling the first and the second layer velocity calculation
recelev=refdata('get','recelev');
fbtime=refdata('get','fbtime');
fbcoord=refdata('get','fbcoord');
shotcoord=refdata('get','shotcoord');
cvpavg=refdata('get','cvpavg');
nshots=refdata('get','nshots');
v1rec = calcvel(fbtime,fbcoord,shotcoord,cvpavg,nshots,recelev);
refdata('set','v1rec',v1rec);
v2rec = calcvel2(fbtime, fbcoord, shotcoord, nshots, cvpavg, recelev);
refdata('set','v2rec', v2rec);
% Update menus
PMTsetmenus;
