function [hcat pcat] = concat_saf_and_regr49(hsub, prof1, head, prof2,opts)
% concat_saf_and_regr49.m
%
% INPUTS:
%       from SAF subset: head and prof2 structures
%       from regr49:     head and prof1 structures
% OUTPUTS:
%       concatenated head and prof structures
%
% pipeline: 
%      1. 

% vers2: 13.Nov.2024 CLH: Sergio created a v2 SAF test prof2iles
%                         use WACCM-X to extend prof2iles to 120 km


addpath /asl/matlib/aslutil            % mktemp
addpath /asl/matlib/h4tools            % rtpwrite
addpath /home/chepplew/gitLib/rtp_prod2/util     % rtp_sub_prof2
addpath /home/chepplew/myLib/matlib/convert_gas_units   % toppmv
addpath /home/chepplew/myLib/matlib/readers             % quick_read_afgl
addpath /home/chepplew/myLib/matlib/sergio_saber_o/

cd /home/chepplew/projects/sarta/matlabcode

% Check inputs: opts
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


% -------------------------------------------------------------
% Append these 49 regression prof2iles to the SAF subset, type 7
% -------------------------------------------------------------
% Need to deal with missing fields between prof1 and prof2
fldnams1 = fieldnames(prof1);
fldnams2 = fieldnames(prof2);
diff12   = setdiff(fldnams1, fldnams2);
diff21   = setdiff(fldnams2, fldnams1);
flddiffs = [diff12; diff21];
% remove unwanted fields diff21 from prof2
prof2 = rmfield(prof2,'pobs');
prof2 = rmfield(prof2,'upwell');
prof2 = rmfield(prof2,'rfreq');

% add diff12 fields to prof2 (49)
prof2.clear  = ones(1,49);
prof2.cc     = zeros(147,49);
prof2.ciwc   = zeros(147,49);
prof2.clwc   = zeros(147,49);
prof2.rcalc  = zeros(1,49);
prof2.robs1  = zeros(1,49);

% add wanted diff21 to prof1 (300,600)
prof1.cfrac   = zeros(size(prof1.rlat));
prof1.salti   = zeros(size(prof1.rlat));
prof1.gas_5   = zeros(size(prof1.gas_1));
prof1.gas_9   = zeros(size(prof1.gas_1));
prof1.scanang = zeros(size(prof1.rlat));
prof1.zobs    = 705000.0*ones(size(prof1.rlat));

% emis, efreq, rho need to be on same frequency grid
emisi = interp1(prof2.efreq(:,1), prof2.emis(:,1), prof1.efreq(:,1),'linear','extrap');
rhoi  = interp1(prof2.efreq(:,1), prof2.rho(:,1), prof1.efreq(:,1),'linear','extrap');
prof2.efreq  = repmat(prof1.efreq(:,1),1,49);
prof2.emis   = repmat(emisi, 1,49);
prof2.rho    = repmat(rhoi, 1,49);

% --------------------------------------------------------------------
% to use rtp_cat_prof2: need field names in same order for prof2 and prof1
prof2 = orderfields(prof2);
prof1 = orderfields(prof1);

% repeat the fieldname checks above,   then concatenate :-
% put the 49 regr first
pcat = rtp_cat_prof(prof1, prof2);
ncat = length(pcat.rlat);

% Update head structure gas field values
%
hcat = headnew;
hcat.glist  = [1,2,3,5,6,9,34]';
hcat.ngas   = 7;
hcat.gunit  = [10,10,10,10,10,10,10]';

% directory to save (if option chosen)
sav.dir = '/home/chepplew/data/sarta/SAF/';
% and a useful comment for the data file:
comment = ['see: /home/chepplew/projects/sarta/matlabcode/' ...
       'subset_saf,load_regr49, splice_WACCM-X, scripts '];

% ---------------------------------------------------------------
% Save mat file 
%
if (mat_save)
  sav.fn = 'saf_sub300_reg49_plv147_400ppm.mat';
%%  sav.fn = 'saf_sub600_reg49_plv147_400ppm.mat';
  if(exist([sav.dir sav.fn]))
     disp('file already exists - not saving')
  else
     disp(['saving data to mat file: ' [sav.dir sav.fn]]);
     save([sav.dir sav.fn],'hcat','pcat','comment','-v7.3');
  end
end
% --------------------------------------------------------------
%     option to save netCDF file for Manuel
if(nc_save)
  mon    = month(tai2dtime(pcat.rtime));
  nprofs = length(pcat.rtime)

  sav.nc = 'saf_sub300_reg49_plv147_400ppm.nc';
%%  sav.nc = 'saf_sub600_reg49_plv147_400ppm.nc;'
  disp(['saving data to netcdf file: ' [sav.dir sav.nc]]);

  ncid      = netcdf.create([sav.dir sav.nc], "NOCLOBBER");
  dimid_np  = netcdf.defDim(ncid,"nprof_dim", nprofs);
  dimid_nl  = netcdf.defDim(ncid,"nlev_dim",147);
  varid_lat = netcdf.defVar(ncid,"latitude","double",dimid_np);
  varid_lon = netcdf.defVar(ncid,"longitude","double",dimid_np);
  varid_mon = netcdf.defVar(ncid,"month","double",dimid_np);
  varid_pre = netcdf.defVar(ncid,"pressure","double", [dimid_nl,dimid_np]);
  varid_alt = netcdf.defVar(ncid,"altitude","double", [dimid_nl,dimid_np]);
  varid_tem = netcdf.defVar(ncid,"temperature","double", [dimid_nl,dimid_np]);
  varid_h2o = netcdf.defVar(ncid,"h2o","double", [dimid_nl,dimid_np]);
  varid_co2 = netcdf.defVar(ncid,"co2","double", [dimid_nl,dimid_np]);
  varid_o3  = netcdf.defVar(ncid,"ozone","double", [dimid_nl,dimid_np]);
  varid_aox = netcdf.defVar(ncid,"atomic_o","double", [dimid_nl,dimid_np]);
  varid_co  = netcdf.defVar(ncid,"co","double", [dimid_nl,dimid_np]);

  varid_str1  = netcdf.getConstant("NC_GLOBAL");
  attname1    = "creation_time";
  attvalue1   = string(datetime("now"));
  varid_str2  = netcdf.getConstant("NC_GLOBAL");
  attname2    = "comments";
  attvalue2   = comment;

  netcdf.endDef(ncid)

  netcdf.putVar(ncid,varid_lat, pcat.rlat)
  netcdf.putVar(ncid,varid_lon, pcat.rlon)
  netcdf.putVar(ncid,varid_mon, mon)
  netcdf.putVar(ncid,varid_pre, pcat.plevs)
  netcdf.putVar(ncid,varid_alt, pcat.palts)
  netcdf.putVar(ncid,varid_tem, pcat.ptemp)
  netcdf.putVar(ncid,varid_h2o, pcat.gas_1)
  netcdf.putVar(ncid,varid_co2, pcat.gas_2)
  netcdf.putVar(ncid,varid_o3,  pcat.gas_3)
  netcdf.putVar(ncid,varid_aox, pcat.gas_34)
  netcdf.putVar(ncid,varid_co,  pcat.gas_5)
  netcdf.reDef(ncid)
  netcdf.putAtt(ncid,varid_str1,attname1,attvalue1)
  netcdf.putAtt(ncid,varid_str2,attname2,attvalue2)

  netcdf.close(ncid);
end


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

% overwrite prof1 fields with new levels data
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
%}

% ============================================================
% EOF
% ============================================================

%{
% --------------------------------------------------------------
% match plevs from 147 levens to the AIRS 101 levels (for ptemp)
%       use first prof2ile as reference.
clear iilvl;
iilvl{1} = find(xprof2new.plevs(:,1) > 0.005  & xprof2new.plevs(:,1) <= 0.1370);
iilvl{2} = find(xprof2new.plevs(:,1) > 0.1370 & xprof2new.plevs(:,1) <= 0.7140);

% get subsets of predictors used in nonLTE
ptemp1 = mean(xprof2new.ptemp(iilvl{1}, irand),1);
ptemp2 = mean(xprof2new.ptemp(iilvl{2}, irand),1);

ptempi1 = mean(prof1.ptemp(1:5,:),1);
ptempi2 = mean(prof1.ptemp(6:9,:),1);

solzen = xprof2new.solzen(irand);
satzen = xprof2new.satzen(irand);
%
% compare layers values to latest regr49.op data used:
%   -> need to trim down number of levels for klayer MAXLEV=120
% ---------------------------------------------------------------
prof101 = prof2;
prof101.ptemp = interp1(plevs147, prof2.ptemp, plevs101, 'linear','extrap');
prof101.gas_1 = interp1(plevs147, prof2.gas_1, plevs101, 'linear','extrap');
prof101.gas_2 = interp1(plevs147, prof2.gas_2, plevs101, 'linear','extrap');
prof101.gas_3 = interp1(plevs147, prof2.gas_3, plevs101, 'linear','extrap');
prof101.gas_5 = interp1(plevs147, prof2.gas_5, plevs101, 'linear','extrap');
prof101.gas_6 = interp1(plevs147, prof2.gas_6, plevs101, 'linear','extrap');
prof101.gas_9 = interp1(plevs147, prof2.gas_9, plevs101, 'linear','extrap');
prof101.nlevs = 98*ones(size(prof2.nlevs));
prof101.plevs = repmat(plevs101,49,1)';
%


rtpwrite(fn.ip, head49,[],prof101,[]);
command = [klayersexe ' fin=' fn.ip ' fout=' fn.op ' > /home/chepplew/logs/klayers/klout.log'];
tic; system(command); toc

[hdo,~,pdo,~]     = rtpread(fn.op);
[hdr49,~,pdr49,~] = rtpread(fn.r49_op);
%}
%{
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
rtpwrite([sav.dir sav.fn], headnew,[],prof1,[]);

sav.fn = 'saf_sub500_plvl147.mat';
prof2_sub = rtp_sub_prof2(xprof2new, irand);
save([sav.dir sav.fn], 'headnew', 'prof2_sub','-v7.3'); 

% ----------------------------------------------------
phome = '/home/chepplew/projects/sarta/figs/saf_subset/';
% cycle: gas_1, gas_2, gas_3, gas_5, gas_6, gas_9, gas_34, ptemp
semilogy(pcat.ptemp, pcat.plevs, '.-','color',[0.6 0.6 0.6])
  grid on; set(gca,'ydir','reverse');xlabel('ppmv');ylabel('hPa');
  title('saf sub300 + regr49 ptemp')
  saveas(gcf,[phome 'saf_sub300_ptemp.png'],'png')

plot(pcat.whichbunch,'.');grid on
  xlabel('prof2 #'); ylabel('atmos group #')
  title('saf sub600 atmos group')
  saveas(gcf, [phome 'saf_sub600_atmos_group.png'],'png');

simplemap(pcat.rlat, pcat.rlon, pcat.whichbunch,2.5)
  title('saf sub300 + regr49 atmos. group')
  saveas(gcf,[phome 'saf_sub300_atmos_group_map.png'],'png')

%}
