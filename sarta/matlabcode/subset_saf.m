function [psub] = subset_saf(opts)
% subset_saf.m
%
% SAF704 new profiles has 21590 profiles on 147 levels. Profiles are ordered:
%    TOA->SFC wherease regr49 are ordered SFC->TOA so reverse regr set.
% Append original Scott 49 regression profiles.
% Splicing based on: add_othergases_arb_pressures.m
%
% pipeline: 
%      1.  Load the SAF LEVELS profiles (NB gunit=10)
%      2.  Subset to 600 or 300
%      2.b repair gas_6 and add to head.glist
%      2.c interpolate cc, ciwc, clwc to the same 147 plevs as gas profiles.
%      3. add SABER atmonic O.
%
% INPUT: opts (options structure) with fields:
%        subset:  [1,2,3,4] option
%        nc_save: save netcdf file [0,1] no/yes
%        mat_save: save matfile [0,1]  no/yes
%
% vers2: 13.Nov.2024 CLH: Sergio created a v2 SAF test profiles
%                         use WACCM-X to extend profiles to 120 km


addpath /asl/matlib/aslutil            % mktemp
addpath /asl/matlib/h4tools            % rtpwrite
addpath /home/chepplew/gitLib/rtp_prod2/util     % rtp_sub_prof
addpath /home/chepplew/myLib/matlib/convert_gas_units   % toppmv
addpath /home/chepplew/myLib/matlib/readers             % quick_read_afgl
addpath /home/chepplew/myLib/matlib/sergio_saber_o/

cd /home/chepplew/projects/sarta/matlabcode

% ==========================================================
% Check input opts structure
if(~exist('opts'))
  disp('Error: need options input')
  return;
end

if(~isfield(opts,'subset'))
  disp('Error: need option subset, 1,2,3 or 4')
  return
end
if(~ismember(opts.subset,[1:4]))
  disp('Error: invalid option for subset')
  return
end
if(~isfield(opts,'nc_save'))
  disp('not saving netcdf data')
  nc_save = false;
else
  if(ismember(opts.nc_save,[0 1])) nc_save = opts.nc_save; end
end
if(~isfield(opts,'mat_save'))
  disp('not saving matlab data')
  mat_save = false;
else
  if(ismember(opts.mat_save,[0 1])) mat_save = opts.mat_save; end
end

% ===============================================================
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

switch opts.subset
  case 1
% subsetting options 1:
% ---------------------
% apply random subset to full set
irand = randsample([1:nprofs], 500);
  case 2
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
  case 3
% subsetting options 3:
% ---------------------
% random over ocean
%   need landfrac for ocean subset

  case 4
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

end     % END switch opts.subset

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
  %saber_O_prof(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk));
  saber_O_prof(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk),'linear','extrap');
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
LPCAT = false;
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

% ----------------------------------
%  fill top part of short profiles
% ----------------------------------
clear prmin iimin
psub.plevs(psub.plevs == 0) = NaN;

% for short profiles extend using afgl (NB profiles are SFC->TOA)
for ip = 1:size(psub.rlat,2)
  [prmin(ip), iimin(ip)] = min(psub.plevs(:,ip));
  if(prmin(ip) > 4.0E-5)
% find best matched AFGL zone to use assigned to IZ
    iipA = find(afgl.a1.pstd <= p49min(ip));
    %%iipB = find(prof49.plevs(:,ip) < p49min(ip),1);
    %%iipB = find(prof49.plevs(:,ip) > p49min(ip),1,'last');
    iipB = i49min(ip)-1;         % ! -1

    [~, IZ] = min( abs(afgl.tstd(iipA(1),:) - prof49.ptemp(i49min(ip),ip) ));
    fprintf(1,'ip, pmin(ip): iipB,IZ %4i, %5.2e, %4i,%4i\n', ...
            ip, p49min(ip), iipB,IZ)


% ----------------------------------
% check for NaN and -ve gas amounts
% ----------------------------------
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
%ineg = find(psub.gas_5 < 0);
%if(~isempty(ineg))
%  psub.gas_5(ineg) = 0.0;
%end
ineg = find(psub.gas_6 < 0);
if(~isempty(ineg))
  psub.gas_6(ineg) = 0.0;
end
%ineg = find(psub.gas_9 < 0);
%if(~isempty(ineg))
%  psub.gas_9(ineg) = 0.0;
%end

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

if(mat_save)
% ---------------------------------------------------------------
%  option to Save mat file 
%
  comment = 'see: /home/chepplew/projects/sarta/matlabcode/subset_saf.m';
  sav.dir = '/home/chepplew/data/sarta/SAF/';
  sav.fn = 'saf_sub300_reg49_plv147_400ppm.mat';
  sav.fn = 'saf_sub600_reg49_plv147_400ppm.mat';
  save([sav.dir sav.fn],'hcat','pcat','comment','-v7.3');
end
% ---------------------------------------------------------------
if(nc_save)
%  option to save netCDF file for Manuel

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

end
% --------------------------------------------------------------
disp('END of subset_saf')

% <<<<<<<<<<<<< EOF >>>>>>>>>>>>>>>>





%{
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
%}
%{
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

rtpwrite(fn.ip, hcat,[],prof101,[]);
command = [klayersexe ' fin=' fn.ip ' fout=' fn.op ' > /home/chepplew/logs/klayers/klout.log'];
tic; system(command); toc

%}
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
