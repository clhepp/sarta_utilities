% compare_rdown_thermal.m

% purpose: compare thermal coefficients and thermal calculations from sarta.

cd /home/chepplew/projects/sarta/cris_hr/

FNATD_long='/asl/s2/hannon/AIRS_prod08/Fit_ftc/Fit_therm/Data/thermdata_m140_long.mat'
FNATD_short='/asl/s2/hannon/AIRS_prod08/Fit_ftc/Fit_therm/Data/thermdata_m140_short.mat';
ARTD_long=load(FNATD_long);
ARTD_short=load(FNATD_short);

fairs=[ARTD_long.freq; ARTD_short.freq];

inname = '/asl/s1/chepplew/projects/sarta/cris_hr/thermdata_cris_hrg4';
CSD=load(inname);
fcris_2235=CSD.freq;   fc=fcris_2235;

FNAD='/asl/s2/hannon/AIRS_prod08/Fit_ftc/Fit_therm/Data/rdown_m130.mat';
FNCD='/asl/s1/chepplew/projects/sarta/cris_hr/rdown_cris_hrg4.mat';
ARD=load(FNAD);
CRD=load(FNCD);

figure(1);clf;plot(fairs,ARD.rdown(:,1),'-',fc,CRD.rdown(:,1),'-');xlim([640 1190]);
  grid on;title('AIRS and CrIS Rdown prof 1 ang 1'); legend('AIRS','CrIS');
  xlabel('wavenumber'); ylabel('radiance');
  % saveas(gcf,'./figs/rdown_airs_cris_prof1_ang1.png','png');



% ---------------------------------------------------------------------
% compare sarta thermal calcs using sarta_crisg4_nov09_wcon_nte (std.res g4) 
% ---------------------------------------------------------------------
addpath /asl/matlib/h4tools
warning 'off'

inFH = fopen('/asl/data/sarta_database/Data_CrIS_apr09/Solar/solardatag4.txt','r');
SD = textscan(inFH,'%f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'CommentStyle','!');
fclose(inFH);

fnrtp='/home/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm.op.rtp';
fnrtp='/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_1013mb_seaemis.rtp';
[head hatt prof patt] = rtpread(fnrtp);
head.nchan = 1329;
head.ichan = int32(SD{1});
head.vchan = single(SD{2});

fortp='/asl/s1/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_1329.op.rtp';
fortp='/asl/s1/chepplew/projects/sarta/cris_hr/regr49_6ang_1013mb_seaemis_1329.op.rtp';
rtpwrite(fortp,head,hatt,prof,patt);


% needs a 1329 channel CrIS file
SARTA_EXE1 = '/asl/packages/sartaV108/BinV201/sarta_crisg4_nov09_wcon_nte';
FIN=/asl/s1/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_1329.op.rtp

