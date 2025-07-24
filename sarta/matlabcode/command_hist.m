% 686 profiles 49 X 7 x 2

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil


srcdr  = '/home/chepplew/data/sarta/prod_2019/generic/';
srcrtp = [srcdr 'r49_1100_98lev_400p_unitemis_seaemis_7angs_night.rtp']; % or _1013_

[head hatt prof patt] = rtpread(srcrtp);
subs.ang = unique(prof.satzen);
idx = struct;
for i=1:7
  idx.a{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) == 1);
  idx.b{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) < 1);
end
subs.names = {'scan angles x 7','2 x surface emissivity: unit and sea'};


kp1 = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
       'REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/RADTEST_SARTA_ChrisH/'];
       
kp2 = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' ...
       'JUNK/Chris_Feb14_2020/'];

kd1 = dir([kp1 'convolved_kcarta_*.mat']);
kd2 = dir([kp2 'individual_prof_convolved_kcarta_crisHI_crisMED*.mat']);

kcr1 = [];
for fn = 1:length(kd1)
  kc1  = load([kd1(fn).folder '/' kd1(fn).name]);
  kcr1 = [kcr1 kc1.med_rcris_all];
  if(fn = 1)
    freq = kc1.med_fcris;
  end
  %fprintf(1,'.')
end

kcr2 = [];
for fn = 1:length(kd2)
  kc2  = load([kd2(fn).folder '/' kd2(fn).name]);
  kcr2 = [kcr2 kc2.med_rcris_all];
  if(fn == 1)
    freq = kc2.med_fcris;
  end
  %fprintf(1,'.')
end




kbt1    = rad2bt(freq, kcr1);
kbt1_mn = nanmean(kbt1,2);
kbt1_sd = nanstd(kbt1,0,2);

kbt2    = rad2bt(freq, kcr2);
kbt2_mn = nanmean(kbt2,2);
kbt2_sd = nanstd(kbt2,0,2);





