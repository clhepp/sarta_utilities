% subset_saf_for_nonlte_modeling.m
%
% SAF704 new profiles has 21590 profiles on 147 levels. Profiles are ordered:
%    TOA->SFC wherease regr49 are ordered SFC->TOA so reverse regr set.
% Append original Scott 49 regression profiles.
% Splicing based on: add_othergases_arb_pressures.m
%
% pipeline: 
%      1.  Load the SAF LEVELS profiles (NB gunit=10)
%      2.  Subset to 600 or 300
%      2.b correct xprofnew.gas_6 and add to headnew.glist
%      2.c interpolate cc, ciwc, clwc to the same 147 plevs as gas profiles.
%      3.  Load the regr49 LEVELS profiles (NB gunit=21)
%      3.a convert mmr to ppmv
%      3.b update CO2ppm and gas_2 amounts to 400 ppmv.
%      3.c  interpolate to the 147 levels of the SAF.
%      4.  Append the regr49 to the SAF subset
%      4.b match the fields of the two sets, adding/removing where different.
%      4.c reorder the fields in preparation for concatenation
%      4.d concatenate the two profile structures.
%      5.  Update head structure with new gas lists
%      6.  save the .mat file
%      7. Interpolate to the airs levels to run klayers to check results
%      8. plot distributions to visualize full field range values.

% vers2: 13.Nov.2024 CLH: Sergio created a v2 SAF test profiles
%                         use WACCM-X to extend profiles to 120 km


addpath /asl/matlib/aslutil            % mktemp
addpath /asl/matlib/h4tools            % rtpwrite
addpath /home/chepplew/gitLib/rtp_prod2/util     % rtp_sub_prof
addpath /home/chepplew/myLib/matlib/convert_gas_units   % toppmv
addpath /home/chepplew/myLib/matlib/readers             % quick_read_afgl
addpath /home/chepplew/myLib/matlib/sergio_saber_o/

cd /home/chepplew/projects/sarta/matlabcode

% original data:
d.home = '/asl/rta/ECMWF_RTTOV_91_25000_and_704_Profiles/';
d.fin  = [d.home 'SAF704_2024/test_profilesV2.mat'];

% Sergio's new version 2 of test_profiles
LV2 = true;

% expect structs: headnew, xprofnew (nprofs=21590, ngas=4)! [1,2,3,6,34 
load(d.fin)
plevs147 = xprofnew.plevs(:,1);

% check contents
fldnams = fieldnames(xprofnew);
nprofs = size(xprofnew.plat,2);

% subsetting options 1:
% ---------------------
% apply random subset to full set
irand = randsample([1:nprofs], 500);

% subsetting options 2:
% ---------------------
% apply random subset to each of the 5 zones
%  try two subsets at 50*5 and 100*5
%    original .whichbunch looks wrong
iiATMo=get_afgl_atmos_type(xprofnew);
irand = [];
nsub  = 50;
for j = 1:6
   %iiwnt = find(xprofnew.whichbunch == j);
   iiwnt = find(iiATMo == j);
   ijunk = randsample([1:length(iiwnt)], nsub);
   irand = [irand iiwnt(ijunk)];
end
nsub = nsub * 6;
psub = rtp_sub_prof(xprofnew, irand);
%%psub.whichbunch = xprofnew.whichbunch(irand);
psub.whichbunch = iiATMo(irand);

% subsetting options 3:
% ---------------------
% random over ocean
%   need landfrac for ocean subset


% subsetting option 4:
% --------------------------------
% random and salti less than 100 m
%   need usgs
addpath /home/chepplew/gitLib/rtp_prod2/util/     % usgs_deg10_dem.m
[salti, landfrac] = usgs_deg10_dem(xprofnew.rlat, xprofnew.rlon);
iiATMo = get_afgl_atmos_type(xprofnew);
irand = [];
nsub  = 50;
for j = 1:6
   iiwnt = find(iiATMo == j & salti < 100);
   ijunk = randsample([1:length(iiwnt)], nsub);
   irand = [irand iiwnt(ijunk)];
end
nsub = nsub * 6;
psub = rtp_sub_prof(xprofnew, irand);
psub.whichbunch = iiATMo(irand);
psub.salti      = salti(irand);
psub.landfrac   = landfrac(irand);

% ---------------------------------------------------------------
% for splicing in missing upper profiles use afgl, load in now
% ---------------------------------------------------------------
afgl.a1 = quick_read_afgl(1,1);
afgl.a2 = quick_read_afgl(1,2);
afgl.a3 = quick_read_afgl(1,3);
afgl.a4 = quick_read_afgl(1,4);
afgl.a5 = quick_read_afgl(1,5);
afgl.a6 = quick_read_afgl(1,6);
afgl.tstd = [afgl.a1.tstd, afgl.a2.tstd, afgl.a3.tstd, afgl.a4.tstd, ...
             afgl.a5.tstd, afgl.a6.tstd];


% ---------------------------------------------------------------
% add gas_34 (atomic O) using SABER data
% ---------------------------------------------------------------
% get afgl atomicO (gas.34)
afgl34 = quick_read_afgl(34,1);

% get Saber atomicO
fn.saber = ['/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/' ...
            'AtomicO_SABER_Mlynczak/atox_athy_night_YY2014_V1.02.nc'];
[atomicO,~]  = read_netcdf_h5(fn.saber);
atomicO.ppmv = toppmv(atomicO.pressure*ones(1,length(atomicO.lat)), ...
               atomicO.ktemp,atomicO.qatox,16,12);
atomicO.lon = wrapTo180(atomicO.lon);

% Match time/lat/lon
daysINyear  = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdaysINyear = cumsum(daysINyear);
for ii = 1 : 12
  boo = find(atomicO.day > cdaysINyear(ii) & atomicO.day <= cdaysINyear(ii+1));
  atomicO.month(boo) = ii;
end

% ------------------------------------------------------
% repair psub gas_34
% ------------------------------------------------------
clear saber_*
mons = month(tai2dtime(psub.rtime));
for ii=1:12
  booO = find(atomicO.month == ii);
  booP = find(mons == ii);
  [A,B,C,D,E] = haversine3_matrix(psub.rlat(booP), psub.rlon(booP), ...
                atomicO.lat(booO), atomicO.lon(booO), 5E7);
  for jj = 1:length(booP)
      [saber_O(:,booP(jj)), saber_p] = add_afgl_g34(atomicO.ppmv(:,booO(B(jj))), ...
            atomicO.pressure, afgl34);
  end
end
for ii = 1 : length(psub.rlat)
  pjunk = psub.plevs(:,ii);
  saber_O_prof(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk));
end
% Rplace NaN with 0.0
clear inan
inan = find(isnan(saber_O_prof));
if(~isempty(inan))
  saber_O_prof(inan) = 0.0;
end
psub.gas_34 = saber_O_prof;

% --------------------------------------------
% repair gas_6 (CH4)
% NB if applying after step 5 use pcat.field(:,50:end) 
% ---------------------------------------------- 
% get afgl profile for region1-6 (TRP,MLS,MLW,SAS,SAW,STD)
LPCAT = true;
if(LPCAT)
  psub = rtp_sub_prof(pcat,[50:299]);
end
iiATM = get_afgl_atmos_type(psub);

% get afgl CH4 (gas_6) and interpolate to plev147
clear gas_6;
for iprof=1:nsub
  afgl6 = quick_read_afgl(6,iiATM(iprof));
  gas_6(:,iprof) = interp1(afgl6.piAtm, afgl6.qiAtm, plevs147, 'linear');
end
% check for NaN and replace with zero
inan = find(isnan(gas_6));
if(~isempty(inan))
  gas_6(inan) = 0.0;
end
% replace psub.gas_6
psub.gas_6 = gas_6;

% --------------------------------------------------
% update co2ppm to 400 ppm (mmr 611.11 g/g dry air)
% --------------------------------------------------
co2_fac = 400.0/max(psub.gas_2(90,:));      % ppmv
max(psub.gas_2(90,:))
psub.co2ppm = 400.0*ones(size(psub.rlat));
psub.gas_2  = co2_fac * psub.gas_2;
% check for NaN and -ve gas amounts
inan = find(isnan(psub.gas_2));
if(~isempty(inan))
  psub.gas_2(inan) = 0.0;
end
ineg = find(psub.gas_1 < 0);
if(~isempty(ineg))
  psub.gas_1(ineg) = 0.0;
end
ineg = find(psub.gas_2 < 0);
if(~isempty(ineg))
  psub.gas_2(ineg) = 0.0;
end
ineg = find(psub.gas_3 < 0);
if(~isempty(ineg))
  psub.gas_3(ineg) = 0.0;
end
ineg = find(psub.gas_5 < 0);
if(~isempty(ineg))
  psub.gas_5(ineg) = 0.0;
end
ineg = find(psub.gas_6 < 0);
if(~isempty(ineg))
  psub.gas_6(ineg) = 0.0;
end
ineg = find(psub.gas_9 < 0);
if(~isempty(ineg))
  psub.gas_9(ineg) = 0.0;
end

% --------------------------------------------------------------
% put cc, ciwc, clwc which are on 137 height grid onto 147 plevs
%     same as other gases. Force upper atmospheric laevels to zero.
x = importdata('/home/chepplew/myLib/data/LP137.txt');
plevs137 = x.data(:,4);
pmin137 = min(plevs137);
imin147 = find(plevs147 < pmin137,1);

cci   = interp1(plevs137, psub.cc, plevs147, 'linear');
ciwci = interp1(plevs137, psub.ciwc, plevs147, 'linear');
clwci = interp1(plevs137, psub.clwc, plevs147, 'linear');
% check for NaN
inan = find(isnan(cci));
if(~isempty(inan))
  cci(inan) = 0.0;
end
inan = find(isnan(ciwci));
if(~isempty(inan))
  ciwci(inan) = 0.0;
end
inan = find(isnan(clwci));
if(~isempty(inan))
  clwci(inan) = 0.0;
end

psub.cc     = cci;
psub.ciwc   = ciwci;
psub.clwc   = clwci;

% --------------------------------------------------------
%  assign psub profiles to nearest AFGL 5-zones (not STD.)
% test ptemp at mesopause (P= 0.00172 to 0.0044 hPa)
Ltesta = find(afgl.a1.pstd >= 0.0017 & afgl.a1.pstd <= 0.0045);
Ltestb = find(plevs147 >= 0.0017 & plevs147 <= 0.0045);

Ta_mn = mean(afgl.tstd(Ltesta,:),1);
Tb_mn = mean(psub.ptemp(Ltestb,:),1);
clear ttmin iimin
for ip = 1:nsub
  [ttmin(ip), iimin(ip)] = min(abs(Tb_mn(ip) - Ta_mn));
end


% --------------------------------------------------------
% add regr49 reference profiles - interpolated to plevs147
% NB:  plevs147 TOA->SFC, prof49.plevs SFC->TOA => reverse 
% --------------------------------------------------------
fn.r49_op = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
           'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];
fn.r49_ip = ['/asl/s1/sergio/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
           'REGR49_H2012_June2016/regr49_1100.ip.rtp'];
% !! check head.gunit :  g/g 21. !!
% NB Only the first 50 levels has values !!
[head49,~,prof49,~] = rtpread(fn.r49_ip);

% load RMM for gases
rmms = importdata('/home/chepplew/myLib/data/RMM_of_HITRAN_gases.txt');

% convert gas units from 10 (mmr) to 21) ppmv
for ii = 1:length(head49.glist)
  gasID = head49.glist(ii);
  gasUI = head49.gunit(ii);
  gasRM = rmms.data(gasID);
  gasAM = prof49.(['gas_' num2str(gasID)]);
  y     = toppmv(prof49.plevs, prof49.ptemp, gasAM, gasRM, gasUI);
  prof49.(['gas_' num2str(gasID)]) = y;
end

% update co2ppm to 400 ppm (mmr 611.11 g/g dry air)
% co2_fac = 6.111e-4/max(prof49.gas_2(:));      % 6.111/5.62   (*10-4  g/g)
co2_fac = 400.0/max(prof49.gas_2(:));      % ppmv
max(prof49.gas_2(:))
prof49.co2ppm = 400.0*ones(size(prof49.co2ppm));
prof49.gas_2  = co2_fac * prof49.gas_2;

% Set plon to 160.0 deg, get salti and landfrac then force to < 100 m
prof49.plon = 160.0*ones(1,49);
[salti, landfrac] = usgs_deg10_dem(prof49.plat, prof49.plon);
iibig = find(salti > 100.0)
if(~isempty(iibig))
  salti(iibig) = 100.0;
end
prof49.salti    = salti;
prof49.landfrac = landfrac;
% set spres = 1013.25 hPa.
prof49.spres    = 1013.25*ones(1,49);

% -------------------------------------------------
% prof49 profiles are mixture of lengths and ranges
%        get min plevs and splice on afgl if needed
% -------------------------------------------------
clear p49min i49min
prof49.plevs(prof49.plevs == 0) = NaN;

% for short profiles extend using afgl (NB profiles are SFC->TOA)
for ip = 1:49
  [p49min(ip), i49min(ip)] = min(prof49.plevs(:,ip));
  if(p49min(ip) > 4.0E-5)
% find best matched AFGL zone to use assigned to IZ
    iipA = find(afglX.pstd <= p49min(ip));
    %%iipB = find(prof49.plevs(:,ip) < p49min(ip),1);
    %%iipB = find(prof49.plevs(:,ip) > p49min(ip),1,'last');
    iipB = i49min(ip)-1;         % ! -1

    [~, IZ] = min( abs(afgl.tstd(iipA(1),:) - prof49.ptemp(i49min(ip),ip) ));
    fprintf(1,'ip, pmin(ip): iipB,IZ %4i, %5.2e, %4i,%4i\n', ...
            ip, p49min(ip), iipB,IZ)

    eval( ['Tjunk=afgl.a' num2str(IZ) '.tstd(:);'] );
    eval( ['Pjunk=afgl.a' num2str(IZ) '.pstd(:);'] );
    T_off = prof49.ptemp(iipB+1,ip) - Tjunk(iipA(1));
    prof49.ptemp(iipB+[1:length(iipA)],ip) = Tjunk(iipA) + T_off;
    prof49.ptemp(prof49.ptemp == 0) = NaN;

    afgl.g1 = quick_read_afgl(1,IZ);
    g1_off  = prof49.gas_1(iipB+1,ip) - afgl.g1.qiAtm(iipA(1));
    prof49.gas_1(iipB+[1:length(iipA)],ip) = afgl.g1.qiAtm(iipA) + g1_off;

    afgl.g2 = quick_read_afgl(2,IZ);
    g2_off  = prof49.gas_2(iipB+1,ip) - afgl.g2.qiAtm(iipA(1));
    prof49.gas_2(iipB+[1:length(iipA)],ip) = afgl.g2.qiAtm(iipA);

    afgl.g3 = quick_read_afgl(3,IZ);
    g3_off  = prof49.gas_3(iipB+1,ip) - afgl.g3.qiAtm(iipA(1));
    prof49.gas_3(iipB+[1:length(iipA)],ip) = afgl.g3.qiAtm(iipA);

    afgl.g5 = quick_read_afgl(5,IZ);
    prof49.gas_5(iipB+[1:length(iipA)],ip) = afgl.g5.qiAtm(iipA);

    afgl.g6 = quick_read_afgl(6,IZ);
    prof49.gas_6(iipB+[1:length(iipA)],ip) = afgl.g6.qiAtm(iipA);

    afgl.g9 = quick_read_afgl(9,IZ);
    prof49.gas_9(iipB+[1:length(iipA)],ip) = afgl.g9.qstd(iipA);
% if the field is extended all other field extended values are set to zero
%    so the test at the start of this loop will fail.
    prof49.plevs(iipB+[1:length(iipA)],ip) = Pjunk(iipA);
    prof49.plevs(prof49.plevs == 0) = NaN;   
    prof49.nlevs(ip)  = length(find(~isnan(prof49.plevs(:,ip))));
  end
  fprintf(1,'.')
end

% --------------------------------------------------
% interpolate to the 147 plevs and revserse TOA->SFC
% --------------------------------------------------
%%plevs147   = plevs147(end:-1:1);
prof       = prof49;
prof.ptemp = []; prof.gas_1 = []; prof.gas_2 = [];
prof.gas_3 = []; prof.gas_5 = []; prof.gas_6 = [];
prof.gas_9 = []; prof.plevs = []; prof.nlevs = [];
prof.plevs = repmat(plevs147,1,49);
for ip = 1:49
  nlev   = prof49.nlevs(ip);
  plevg  = log(prof49.plevs(nlev:-1:1,ip));
  psurf  = prof49.plevs(1,ip);
  prof.ptemp(:,ip) = interp1(plevg, prof49.ptemp(nlev:-1:1,ip), log(plevs147), ...
        [],'extrap');
  prof.gas_1(:,ip) = interp1(plevg, prof49.gas_1(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_2(:,ip) = interp1(plevg, prof49.gas_2(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_3(:,ip) = interp1(plevg, prof49.gas_3(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_5(:,ip) = interp1(plevg, prof49.gas_5(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_6(:,ip) = interp1(plevg, prof49.gas_6(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
  prof.gas_9(:,ip) = interp1(plevg, prof49.gas_9(1:nlev,ip), log(plevs147), ...
         'linear','extrap');
% check range limits. SFC pressure
  prof.gas_1(prof.gas_1 < 0) = 0;
  prof.gas_2(prof.gas_2 < 0) = 0;
  prof.gas_3(prof.gas_3 < 0) = 0;
  iiovr = find(plevs147 > prof49.plevs(1,ip));
  iirng = find(plevs147 <= prof49.plevs(1,ip));
  prof.nlevs(ip)   = length(iirng);
  %prof.plevs(:,ip) = plevs147(1:length(iirng));
  prof.plevs(iiovr,ip) = 0.0;
end

% ----------- assign needed fields to prof (49) ------
prof.rlat   = prof.plat;
prof.rlon   = prof.plon;
prof.rtime  = repmat(dtime2tai('2020/06/01'), 1,49);
prof.whichbunch = 7.0 * ones(1,49);

% ------------------------------------------
% add atomicO to prof49, are all in June
% ------------------------------------------
clear saber_*
mons = month(tai2dtime(prof.rtime));
ii = 6;
  booO = find(atomicO.month == ii);
  booP = find(mons == ii);
  [A,B,C,D,E] = haversine3_matrix(prof.rlat(booP), prof.rlon(booP), ...
                atomicO.lat(booO), atomicO.lon(booO), 5E7);
  for jj = 1:length(booP)
     [saber_O(:,booP(jj)), saber_p] = add_afgl_g34(atomicO.ppmv(:,booO(B(jj))), ...
           atomicO.pressure, afgl34);
  end
% interpolate to prof.plevs
for ii = 1 : length(prof.rlat)
  pjunk = prof.plevs(:,ii);
  saber_O_prof(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk));
end
% replace existing zero values in prof.gas_34
prof.gas_34 = saber_O_prof;

% -------------------------------------------------------------
% Append these 49 regression profiles to the SAF subset, type 7
% -------------------------------------------------------------
% Need to deal with missing fields between psub and prof
fldnams1 = fieldnames(psub);
fldnams2 = fieldnames(prof);
diff12   = setdiff(fldnams1, fldnams2);
diff21   = setdiff(fldnams2, fldnams1);
flddiffs = [diff12; diff21];
% remove unwanted fields diff21 from prof
prof = rmfield(prof,'pobs');
prof = rmfield(prof,'upwell');
prof = rmfield(prof,'rfreq');

% add diff12 fields to prof (49)
prof.clear  = ones(1,49);
prof.cc     = zeros(147,49);
prof.ciwc   = zeros(147,49);
prof.clwc   = zeros(147,49);
prof.rcalc  = zeros(1,49);
prof.robs1  = zeros(1,49);

% add wanted diff21 to psub (300,600)
psub.cfrac   = zeros(size(psub.rlat));
psub.salti   = zeros(size(psub.rlat));
psub.gas_5   = zeros(size(psub.gas_1));
psub.gas_9   = zeros(size(psub.gas_1));
psub.scanang = zeros(size(psub.rlat));
psub.zobs    = 705000.0*ones(size(psub.rlat));

% emis, efreq, rho need to be on same frequency grid
emisi = interp1(prof.efreq(:,1), prof.emis(:,1), psub.efreq(:,1),'linear','extrap');
rhoi  = interp1(prof.efreq(:,1), prof.rho(:,1), psub.efreq(:,1),'linear','extrap');
prof.efreq  = repmat(psub.efreq(:,1),1,49);
prof.emis   = repmat(emisi, 1,49);
prof.rho    = repmat(rhoi, 1,49);

% --------------------------------------------------------------------
% to use rtp_cat_prof: need field names in same order for prof and psub
prof = orderfields(prof);
psub = orderfields(psub);

% repeat the fieldname checks above,   then concatenate :-

pcat = rtp_cat_prof(prof, psub);
ncat = length(pcat.rlat);

% Update head structure gas field values
%
hcat = headnew;
hcat.glist  = [1,2,3,5,6,9,34]';
hcat.ngas   = 7;
hcat.gunit  = [10,10,10,10,10,10,10]';

% ---------------------------------------------------------------
% Save mat file 
%
 comment = 'see: /home/chepplew/projects/sarta/matlabcode/subset_saf_for_nonlte_modeling.m';
 sav.dir = '/home/chepplew/data/sarta/SAF/';
 sav.fn = 'saf_sub300_reg49_plv147_400ppm.mat';
 sav.fn = 'saf_sub600_reg49_plv147_400ppm.mat';
 save([sav.dir sav.fn],'hcat','pcat','comment','-v7.3');

% --------------------------------------------------------------



% --------------------------------------------------------------
%  run through klayers
%    interpolate from 147 to 101 plevs (maxlev=120)
% --------------------------------------------------------------
klayersexe = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
temp.dir = mktemp();
fn.ip    = [temp.dir '_ip.rtp'];
fn.op    = [temp.dir '_op.rtp'];

% trim/interp onto MAXLEV 120 vertical grid for klayers
% load AIRS plevs
x=importdata('/home/chepplew/myLib/data/airs_plevs.txt')
x2 = reshape(x',[],105);
plevs101 = x2(~isnan(x2));
% BUT MUST trnucate plevs101 to 1013mb SFC for interpolation to plevs147 etc
plevs101 = plevs101(1:end-3);
clear x x2;

prof101 = pcat;

ptempi = interp1(pcat.plevs(:,1),pcat.ptemp, plevs101, 'linear');
gas1i  = interp1(pcat.plevs(:,1),pcat.gas_1, plevs101, 'linear');
gas2i  = interp1(pcat.plevs(:,1),pcat.gas_2, plevs101, 'linear');
gas3i  = interp1(pcat.plevs(:,1),pcat.gas_3, plevs101, 'linear');
gas5i  = interp1(pcat.plevs(:,1),pcat.gas_5, plevs101, 'linear');
gas6i  = interp1(pcat.plevs(:,1),pcat.gas_6, plevs101, 'linear');
gas9i  = interp1(pcat.plevs(:,1),pcat.gas_9, plevs101, 'linear');
gas34i = interp1(pcat.plevs(:,1),pcat.gas_34, plevs101, 'linear');
 
cci   = interp1(pcat.plevs(:,1), pcat.cc,   plevs101, 'linear');
ciwci = interp1(pcat.plevs(:,1), pcat.ciwc, plevs101, 'linear');
clwci = interp1(pcat.plevs(:,1), pcat.clwc, plevs101, 'linear');

% overwrite psub fields with new levels data
prof101.nlevs  = 101*ones(size(pcat.nlevs));
prof101.plevs  = repmat(plevs101, ncat,1)';
prof101.ptemp  = ptempi;
prof101.gas_1  = gas1i;
prof101.gas_2  = gas2i;
prof101.gas_3  = gas3i;
prof101.gas_5  = gas5i;
prof101.gas_6  = gas6i;
prof101.gas_9  = gas9i;
prof101.gas_34 = gas34i;
prof101.cc     = cci;
prof101.ciwc   = ciwci;
prof101.clwc   = clwci;

%


rtpwrite(fn.ip, hcat,[],prof101,[]);
command = [klayersexe ' fin=' fn.ip ' fout=' fn.op ' > /home/chepplew/logs/klayers/klout.log'];
tic; system(command); toc


% ============================================================
% EOF
% ============================================================

%{

% --------------------------------------------------------------
% match plevs from 147 levens to the AIRS 101 levels (for ptemp)
%       use first profile as reference.
clear iilvl;
iilvl{1} = find(xprofnew.plevs(:,1) > 0.005  & xprofnew.plevs(:,1) <= 0.1370);
iilvl{2} = find(xprofnew.plevs(:,1) > 0.1370 & xprofnew.plevs(:,1) <= 0.7140);

% get subsets of predictors used in nonLTE
ptemp1 = mean(xprofnew.ptemp(iilvl{1}, irand),1);
ptemp2 = mean(xprofnew.ptemp(iilvl{2}, irand),1);

ptempi1 = mean(psub.ptemp(1:5,:),1);
ptempi2 = mean(psub.ptemp(6:9,:),1);

solzen = xprofnew.solzen(irand);
satzen = xprofnew.satzen(irand);
%
% compare layers values to latest regr49.op data used:
%   -> need to trim down number of levels for klayer MAXLEV=120
% ---------------------------------------------------------------
prof101 = prof;
prof101.ptemp = interp1(plevs147, prof.ptemp, plevs101, 'linear','extrap');
prof101.gas_1 = interp1(plevs147, prof.gas_1, plevs101, 'linear','extrap');
prof101.gas_2 = interp1(plevs147, prof.gas_2, plevs101, 'linear','extrap');
prof101.gas_3 = interp1(plevs147, prof.gas_3, plevs101, 'linear','extrap');
prof101.gas_5 = interp1(plevs147, prof.gas_5, plevs101, 'linear','extrap');
prof101.gas_6 = interp1(plevs147, prof.gas_6, plevs101, 'linear','extrap');
prof101.gas_9 = interp1(plevs147, prof.gas_9, plevs101, 'linear','extrap');
prof101.nlevs = 98*ones(size(prof.nlevs));
prof101.plevs = repmat(plevs101,49,1)';
%


rtpwrite(fn.ip, head49,[],prof101,[]);
command = [klayersexe ' fin=' fn.ip ' fout=' fn.op ' > /home/chepplew/logs/klayers/klout.log'];
tic; system(command); toc

[hdo,~,pdo,~]     = rtpread(fn.op);
[hdr49,~,pdr49,~] = rtpread(fn.r49_op);



% ---------------------------------------------------
% plot distros:
% ---------------------------------------------------
histogram(ptemp1,[170:2:290],'Normalization','probability')
  hold on
  histogram(ptemp2,[170:2:290],'Normalization','probability')
  grid on;  legend('ptemp1','ptemp2')
  title('SAF.21590 subset 500 ptemp1,2')
  %saveas(gcf, './figs/saf_subset_500_ptemp_pdf.fig','fig')

% save the 101 levels rtp
sav.dir = '/home/chepplew/data/sarta/SAF/';
sav.fn  = 'saf_sub500_plv101.rtp';
rtpwrite([sav.dir sav.fn], headnew,[],psub,[]);

sav.fn = 'saf_sub500_plvl147.mat';
prof_sub = rtp_sub_prof(xprofnew, irand);
save([sav.dir sav.fn], 'headnew', 'prof_sub','-v7.3'); 

% ----------------------------------------------------
phome = '/home/chepplew/projects/sarta/figs/saf_subset/';
% cycle: gas_1, gas_2, gas_3, gas_5, gas_6, gas_9, gas_34, ptemp
semilogy(pcat.ptemp, pcat.plevs, '.-','color',[0.6 0.6 0.6])
  grid on; set(gca,'ydir','reverse');xlabel('ppmv');ylabel('hPa');
  title('saf sub300 + regr49 ptemp')
  saveas(gcf,[phome 'saf_sub300_ptemp.png'],'png')

plot(pcat.whichbunch,'.');grid on
  xlabel('prof #'); ylabel('atmos group #')
  title('saf sub600 atmos group')
  saveas(gcf, [phome 'saf_sub600_atmos_group.png'],'png');

simplemap(pcat.rlat, pcat.rlon, pcat.whichbunch,2.5)
  title('saf sub300 + regr49 atmos. group')
  saveas(gcf,[phome 'saf_sub300_atmos_group_map.png'],'png')
  

%}

%{
% ==========================================================
%     option to save netCDF file for Manuel
% ==========================================================

mon    = month(tai2dtime(pcat.rtime));
nprofs = length(pcat.rtime)

sav.nc = 'saf_sub300_reg49_plv147_400ppm.nc';
sav.nc = 'saf_sub600_reg49_plv147.nc';

ncid      = netcdf.create([sav.dir sav.nc], "NOCLOBBER");
dimid_np  = netcdf.defDim(ncid,"nprof_dim", nprofs);
dimid_nl  = netcdf.defDim(ncid,"nlev_dim",147);
varid_lat = netcdf.defVar(ncid,"latitude","double",dimid_np);
varid_lon = netcdf.defVar(ncid,"longitude","double",dimid_np);
varid_mon = netcdf.defVar(ncid,"month","double",dimid_np);
varid_pre = netcdf.defVar(ncid,"pressure","double", [dimid_nl,dimid_np]);
varid_tem = netcdf.defVar(ncid,"temperature","double", [dimid_nl,dimid_np]);
varid_co2 = netcdf.defVar(ncid,"co2","double", [dimid_nl,dimid_np]);
varid_o3  = netcdf.defVar(ncid,"ozone","double", [dimid_nl,dimid_np]);
varid_aox = netcdf.defVar(ncid,"atomic_o","double", [dimid_nl,dimid_np]);
varid_co  = netcdf.defVar(ncid,"co","double", [dimid_nl,dimid_np]);

netcdf.endDef(ncid)

netcdf.putVar(ncid,varid_lat, pcat.rlat)
netcdf.putVar(ncid,varid_lon, pcat.rlon)
netcdf.putVar(ncid,varid_mon, mon)
netcdf.putVar(ncid,varid_pre, pcat.plevs)
netcdf.putVar(ncid,varid_tem, pcat.ptemp)
netcdf.putVar(ncid,varid_co2, pcat.gas_2)
netcdf.putVar(ncid,varid_o3,  pcat.gas_3)
netcdf.putVar(ncid,varid_aox, pcat.gas_34)
netcdf.putVar(ncid,varid_co,  pcat.gas_5)

netcdf.close(ncid);

%}
