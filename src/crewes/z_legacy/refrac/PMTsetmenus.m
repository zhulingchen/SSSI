% This diagram shows the dependancies between 
% process -> data -> menu enabling.
%
% Process	Data output	Menu		Submenu
% -------	-----------	----		-------
% load		shotcoord	Process		rectimecheck
%	 	fbcoord				autopick
% 		fbtimes		
% 		nshots		Display		refracted arrivals
% 		nrecs
% ------------------------------------------------------------------
% rec t. check	diffmat		Display		reciprocal time
% ------------------------------------------------------------------
% autopick	cvpi		Process		autoreject CVP
% 		cvpj				average CVP
%                                               Repick CVP
% 
% 				Display		CVP statistics
% 						CVP average
%                                               CVP Polyfit
% 
% 				Edit		CVPs
% ------------------------------------------------------------------
% averageCVP	cvpavg		Process		calc vel
% 						time analysis
% 
% 				Edit		CVP average
% ------------------------------------------------------------------
% calcvel 	v1rec		Display		velocity model
% 		v2rec		
% 				Edit		velocity model
% 
% 				File		Export - velocity
% ------------------------------------------------------------------
% time analysis	td1		Display		plus time
% 		plust
% 				Process		delay time analysis
%						time rejection
%						depth calc
% ------------------------------------------------------------------
% time rejection plust		
% ------------------------------------------------------------------
% depth calc	depth		process		Static
% 		plustime	
% 
% 				Display		depth model
% 
% 				Edit		depth model
% 
% 				File		Export - depth
% ------------------------------------------------------------------
% static	shotstat	Display		static
% 		recstat		
% 				File		Export - static
function PMTsetmenus
% Get all the menu handles
autorejectm = refdata('get', 'autorejectm');
avgcvpm = refdata('get', 'avgcvpm');
calcvelm = refdata('get', 'calcvelm');
timeanalm = refdata('get', 'timeanalm');
delaytm = refdata('get','delaytm');
avgdepthm = refdata('get', 'avgdepthm');
staticm = refdata('get', 'staticm');
rectimem = refdata('get', 'rectimem');
timerejectm = refdata('get', 'timerejectm');
autopickm = refdata('get', 'autopickm');
dispshotsm = refdata('get', 'dispshotsm');
disprectimem = refdata('get', 'disprectimem');
dispcvpstatm = refdata('get', 'dispcvpstatm');
%dispcvpstat2m = refdata('get', 'dispcvpstat2m');
dispcvpavgm = refdata('get', 'dispcvpavgm');
dispvelm = refdata('get', 'dispvelm');
dispdepthm = refdata('get', 'dispdepthm');
dispstaticm = refdata('get', 'dispstaticm');
editcvpm = refdata('get', 'editcvpm');
editcvpavgm = refdata('get', 'editcvpavgm');
editvelm = refdata('get', 'editvelm');
editdepthm = refdata('get', 'editdepthm');
dispplustm = refdata('get', 'dispplustm');
expvelm = refdata('get', 'expvelm');
expdepthm = refdata('get', 'expdepthm');
expstatm = refdata('get', 'expstatm');
autorejectm = refdata('get', 'autorejectm');
disppolym = refdata('get','disppolym');
repickcvpm = refdata('get','repickcvpm');
setmenustatus([rectimem autopickm dispshotsm], 'shotcoord');
setmenustatus([disprectimem], 'diffmat');
setmenustatus([autorejectm avgcvpm repickcvpm dispcvpstatm dispcvpavgm disppolym editcvpm], 'cvpi');
setmenustatus([calcvelm timeanalm editcvpavgm], 'cvpavg');
setmenustatus([dispvelm editvelm expvelm], 'v1rec');
setmenustatus([timerejectm avgdepthm], 'td1');
setmenustatus([delaytm dispplustm avgdepthm], 'plust');
setmenustatus([staticm dispdepthm editdepthm expdepthm], 'depth');
setmenustatus([dispstaticm expstatm], 'shotstat');
