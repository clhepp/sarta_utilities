% concatenate group of RTPs used by Sergio for 6 sat * 19 sol * 48 prfs
% for extended nonLTE calcs.

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools


rtp.home = ['/home/sergio/KCARTA/NONLTE_PRODUCTION/' ...
   'VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Apr2021_NewNLTEProfiles_SAVETHISGOOD/' ...
   'RTP_100/RTP_100_19solarangs_MLPuertas_Mar2020/'];
   
rtp.list = dir([rtp.home 'rtp_regress_*_400.op.rtp']);

i = 1;
[head,hatt,prof,patt] = rtpread([rtp.list(i).folder '/' rtp.list(i).name]);

for i=2:length(rtp.list)
  [hd,~,pd,~] = rtpread([rtp.list(i).folder '/' rtp.list(i).name]);
  [hda, pda]  = cat_rtp(head, prof, hd, pd);
  
  head = hda;
  prof = pda;
  
  fprintf(1,'%d\t',i)
end

drout = '/home/chepplew/data/sarta/prod_2021/generic/';
fnout = 'r49_400p_6satzen_19solzen_for_nonlte.rtp';

rtpwrite([drout fnout], head, hatt, prof, patt);



