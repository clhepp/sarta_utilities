function [wac] = subset_wacmmx_to_saf(prof,opts)
% splice_waccm_to_saf_r49
%
% load WACCM, subset sample, splce to SAF at ECMWF level
% ECMWF Pmin = 0.02 hPa Desired upper limit = 2.0E-5 hPa
%
%  INPUT:
%       prof structure from subset_saf()
%       opts: structure with fields:
%             mat_save, nc_save
% Notes: WACCM-X profiles: co2vmr, ch4vmr,n20vmr,f11vmr,f12vmr,
%   H, NO, NOX, O, O2, O3, OH, TElec, T, TIon, TS, U, V
%
%
%

addpath /home/chepplew/myLib/matlib/readers

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

% required sample size (# SAF or regr49)
nsam = size(prof.rlat,2);

% splicing pressure (top of ECMWF ~ 0.02 hPa)
p_splice = 0.02;       % hPa
p_toa    = 2.0E-5;     % hPa

% load WACCM and subset to desired sample size
wac.dir = dir(['/home/chepplew/data/WACCM-X/*.nc']);
ifn = 1;
wac.fn  = [wac.dir(ifn).folder '/' wac.dir(ifn).name];

[data att] = read_netcdf_h5(wac.fn);
% [1440x721x273] [ lon x lat x level]

% subset WACCM
%   would be nice to match locations with the prof.rlat/rlon
nwac  = 1440*721;
irand = randsample([1:nwac], nsam);

% latlon = meshgrid(-90:0.25:90, 0:0.25:359.75)
latlon = meshgrid(1:721, 1:1440);
[ilon ilat] = ind2sub(size(latlon),irand);   % irow<=lon, icol<=lat
wac.slat = data.lat([ilat]);
wac.slon = data.lon([ilon]);

% cut altitude range required (ECMWF to TOA)
isel = find(data.lev <= p_splice & data.lev > p_toa);

% collect profiles to pass on: (Gas units VMR dimensionless)
wac.co2vmr = data.co2vmr;        % single value
wac.plev   = data.lev(isel);
  junk = permute(data.T,[3 1 2]);
wac.ptemp  = junk(isel,irand);
  junk = permute(data.O,[3 1 2]);
wac.O      = junk(isel,irand);
  junk = permute(data.O2,[3 1 2]);
wac.O2     = junk(isel, irand);
  junk = permute(data.O3, [3 1 2]);
wac.O3     = 1.0E6*junk(isel, irand);

% splice onto the prof fields: interpolate waccm to prof levels
%     and offset at splice level and scaled (TBD)
for ip = 1:nsam
  iiplev = find(prof.plevs(:,ip) < p_splice,1,'last');
  iiplev = find(prof.plevs(:,ip) > p_toa & prof.plevs(:,ip) <= p_splice);
  tempo = prof.ptemp(iiplev(end),ip);
  tempa = wac.ptemp(end,ip);
  Dtemp = tempo - tempa;
  ijunk = interp1(wac.plev, wac.ptemp(:,ip), prof.plevs(iiplev,ip),...
          'linear','extrap');
  prof.ptemp(iiplev,ip) = ijunk + Dtemp;
% Ozone !! check units !!
  ozo   = prof.gas_3(iiplev(end),ip);
  oza   = wac.O3(end,ip);
  Doz   = ozo - oza;
  ijunk = interp1(wac.plev, wac.O3(:,ip), prof.plevs(iiplev,ip),...
          'linear','extrap');
  prof.gas_3(iiplev,ip) = ijunk + Doz;

end

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

%{
% ==============================================================
% get some analysis - eg PCA
% ==============================================================
ptemp_mn = mean(wac.ptemp,2,'omitnan');
[coeff,score,latent,tsquared,explained,mu] = pca(wac.ptemp./ptemp_mn);
plot(cumsum(latent./sum(latent)),'.-')
grid on; xlim([0 50]);xlabel('component #');ylabel('cummulative contribution');
title('WACCM-X subset 349 ptemp eigenvector');
% saved to: /home/chepplew/projects/WACCM/Figs/waccmx_sub349_ptemp_eigenv

ptemp_mn = mean(prof.ptemp,2,'omitnan');
[coeff,score,latent,tsquared,explained,mu] = pca(prof.ptemp./ptemp_mn);
% using original psub from SAF/AFGL subset:
ptemp_mn = mean(psub.ptemp,2,'omitnan');
[coeff,score,latent,tsquared,explained,mu] = pca(psub.ptemp./ptemp_mn);

%}
