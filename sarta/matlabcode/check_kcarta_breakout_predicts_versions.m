function [] = check_kcarta_breakout_predicts_versions()

% Checking versions of L2S OD breakout predicts:

all_brkouts = {'wvbandFMW','wvbandFM','wvbandF','so2bandS','so2bandF', ...
   'nh3bandNH3','n2ohno3bandN2O','n2ohno3bandHNO3','n2ohno3bandF', ...
   'FO','FW','FWO','FWO_dHDO','FWOP','FWOP_Orig','FWOP_CO2_1.03', ...
   'FWOP_CO2_1.05', 'F', 'FD',' FDO', 'ch4bandCH4', ...
   'ch4bandF','cobandF','cobandFC','cobandFCO','cobandFCOW','cobandFCOWP'};
length(all_brkouts)

dp0='~sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Dec2018_AIRS2834/';
dp1='~sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/';

%dn0=[dp0 'FWOP_Orig/'];
%dn1=[dp0 'FWOP_CO2_1.03/'];
%dn2=[dp0 'FWOP_CO2_1.05/'];
%dn4=[dp1 'FWOP/'];

for ibrk = 21:27;
  this_brkout = all_brkouts{ibrk}

  dn0=[dp0 this_brkout '/'];
  dn1=[dp1 this_brkout '/'];

  fn0=dir([dn0 'convolved_kcarta_*.mat']);
  fn1=dir([dn1 'convolved_kcarta_*.mat']);
  %fn2=[dn2 'convolved_kcarta_FWOP_1.mat'];
  %fn4=[dn4 'convolved_kcarta_FWOP_1.mat'];
  disp(['fn0: ' num2str(length(fn0)) ' fn1: ' num2str(length(fn1))])

  iprf = 1;
  x0(iprf) = load([fn0(iprf).folder '/' fn0(iprf).name]);
  x1(iprf) = load([fn1(iprf).folder '/' fn1(iprf).name]);

  %figure(1);clf;plot(x0.fairs, x0.rairs_all(1,:,20),'-');
  figure(1);clf;
  for j = 1:4:14
   cla; 
   plot(x0.fairs, x0.rairs_all(j,:,20) - x1.rairs_all(j,:,20),'-', ...
        x0.fiasi, x0.riasi_all(j,:,20) - x1.riasi_all(j,:,20),'-', ...
        x0.fcris, x0.rcris_all(j,:,20) - x1.hi_rcris_all(j,:,20),'-');
   title([this_brkout  ' prof: ' num2str(iprf) ' ang: ' num2str(j) ' levl: ' num2str(20)])
   pause(2)
  end
%
end


