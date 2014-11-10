function avgcvpcb
% Function calling the averaging cross over point function
cvpi=refdata('get','cvpi');
cvpj=refdata('get','cvpj');
nshots=refdata('get','nshots');
[cvpavg cvpstd cvpfold]=avgcvp(cvpi,cvpj,nshots);
refdata('set','cvpavg',cvpavg);
% Update menus
PMTsetmenus;
