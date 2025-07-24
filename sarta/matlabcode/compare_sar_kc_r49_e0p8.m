cd /home/chepplew/projects/sarta/cris_hr
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

iang = [1:49; 50:98; 99:147; 148:196; 197:245; 246:294];

dpk=['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/'...
       'JUNK/xconvolved_kcarta_crisHI_xconvolved_kcarta_crisHI_1013mb_0p8emiss.mat'];
dpk=['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/'...
      'JUNK/xconvolved_kcarta_crisHI_xconvolved_kcarta_crisHI_1013mb_0p8emiss_acos35.mat'];      
dpk=['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/'...
      'JUNK/xconvolved_kcarta_crisHI_xconvolved_kcarta_crisHI_1013mb_0p8emiss_integrate.mat'];
KR = load(dpk);
bkc = real(rad2bt(KR.fcris,KR.rcris_all));
for ia = 1:6
  bkcm(:,ia) = nanmean(bkc(:,iang(ia,:)),2);
end

% get the rtp file - and change head structure to run w/sarta CrIS hiRes
fnrtp = 'regr_rtp_6angs_49profs_1013mb_0p8emis.rtp';
load('rtp_head_str_2235g4.mat');                     % <- hdr2 w/CrIS hrg4 channels
[hdr har pdr par] = rtpread(fnrtp);
fortp='regr49_6angs_1013mb_0p8emis_2235g4.rtp';
rtpwrite(fortp, hdr2, har, pdr, par);

%{
% now run sarta   
SARTA_EXE=/home/chepplew/gitLib/osarta/bin/sarta_g4_cris_hrg4_400p_full_th2
SARTA_EXE=/home/chepplew/gitLib/sarta/bin/crisg4_oct16
FIN=/home/chepplew/projects/sarta/cris_hr/regr49_6angs_1013mb_0p8emis_2235g4.rtp
FOUT=/asl/s1/chepplew/projects/sarta/cris_hr/sar_r49_6a_0p8.rtp
%}

fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_r49_6a_0p8.rtp';
fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_r49_6a_sea.rtp';
[hds has pds pas] = rtpread(fnsar);

bsc = real(rad2bt(hds.vchan,pds.rcalc));
for ia = 1:6
  bscm(:,ia) = nanmean(bsc(:,iang(ia,:)),2);
end

fc=hds.vchan;    
figure(1);clf;hold on;grid on;
  for ia=1:6 plot(fc,bkcm(:,ia)-bscm(:,ia),'-'); end
  xlim([640 2550]); %axis([640 1110 -10.0 0.2]); 
title('sar crisg4.oct16 v kc r49 6angs e0p8 int');xlabel('wn cm^{-1}');ylabel('BT bias (K)'); 
 %saveas(gcf,'./figs/sar_crisg4_oct16_r49_kc_bias_6angs_0p8_integrtd_LW.png','png'); 

