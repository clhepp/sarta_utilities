% ------------------------------------
% Compare SARTA TOA radiance w/ kcarta
% ------------------------------------
cd /home/chepplew/projects/sarta/cris_hr/
cd /home/chepplew/projects/sarta/prod_2018/

addpath /asl/matlib/h4tools
addpath /asl/packages/airs_decon/source          % hamm_app.m
addpath /asl/matlib/aslutil                      % rad2bt.m
addpath /asl/matlib/plotutils                    % aslprint
addpath /home/chepplew/projects/sarta/matlabcode
warning 'off'

% Choose which regression set to compare.
all_cregr = {'r49','saf704'};
cregr     = 'saf704';

% Enter which sensor to compare w/kcarta {'AIRS' or 'CRIS'}
csens = 'CRIS';

% Default scan angles (LW) for prod_2018 use first 8.
scang = [0   28.6101   38.5184   45.2214   46.3539   48.3250];
scang = [0    8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
        65.0428   69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];
scang = scang(1:8);
iang  =  [1:numel(scang)];

%  1:              BASELINE TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------
% Get the KLAYERS RTP file used to generate the SARTA and kcarta data
% -------------------------------------------------------------------
kpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018'...
       '/RAD1100_unitemis/'];
if(strcmp(cregr,'r49'))
  fnlst = dir([kpath 'convolved_kcarta_RAD1100_*_radiances.mat']);
% reorder these in profile number order
  for i=1:49 fnparts{i} = strsplit(fnlst(i).name, '_'); 
    prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
  end
  [B IB ] = sort(prfnums);  

  fnrtp = [kpath 'regr49_1100_400ppm.op.rtp';
  fnrtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_400ppm_2235g4.op.rtp';
end
if(strcmp(cregr,'saf704'))
  %/home/sergio/MATLABCODE/REGR_PROFILES/ECMWF_SAF_137Profiles/save_SAF_704_profiles_29-Apr-2016_1100mb.op.rtp
  fnrtp = '/home/chepplew/gitLib/o_sarta/test/SAF704_profs_29Apr2016_2235_1100mb.op.rtp';
  fnrtp = ['/home/chepplew/projects/sarta/cris_hr/SAF704_regrprofs_25May2016_xmb_2235g4.op.rtp'];
end
[hdr har pdr par] = rtpread(fnrtp);

% -------------------------------------------------------------------	   
% get the SARTA calculated radiances from the regression profiles
% -------------------------------------------------------------------
if(strcmp(cregr,'r49'))
  sarfile = '/home/chepplew/gitLib/sarta/test/regr49_hrg4_sar_optr.rtp';
  sarfile = '/home/chepplew/gitLib/sarta/test/regr49_hrg4_sar_optr_MW.rtp';
  sarfile = '/home/chepplew/projects/sarta/cris_hr/regr49_hrg4_sar_wcon_allSets.rtp';
  sarfile = '/home/chepplew/projects/sarta/cris_hr/regr49_hrg4_sar_wcon_wn2o.rtp';
  sarfile = '/home/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_sar_7set.rtp';
  sarfile = '/asl/s1/chepplew/projects/sarta/airs/sar_r49_apr08_m140_wcon_nte.rtp';
  sarfile = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_oct16_r49_1013_420_6a.rtp';
  sarfile = '/asl/s1/chepplew/data/sarta_wrk/cris_hr/sar_saf704_r49_20170202a.rtp';
end
if(strcmp(cregr,'saf704'))
  sarfile = '/home/chepplew/data/sarta_wrk/cris_hr/sar_saf704.rtp';
  sarfile = '/home/chepplew/projects/sarta/cris_hr/sar704_xmb_2235g4.rtp';
end
[fpath fname fext] = fileparts(sarfile);
[hds has pds pas]  = rtpread(sarfile);
 
% ----------------------------------------------------
% get the kcarta reference set for the SAME regr49 set
% ----------------------------------------------------
if(strcmp(cregr,'r49'))
  dpk   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49/RAD1013/';
  dpk   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/RAD1100/';
  dpk   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/RAD1013_CO2x1.05/';
   flist = dir(strcat(dpk,'convolved_kcarta*radiances.mat'));
   for i=1:49 names{i}=flist(i).name; end
   fnum = 49;
   %ropro = [10:19 1 20:29 2 30:39 3 40:49 4:9];   % the order the mat files are loaded.
end
clear names X ndx;
if(strcmp(cregr,'saf704'))
   dpk   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/RAD1100/';
   dpk   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/RAD1013/';
   dpk   = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/'; %
   flist = dir(strcat(dpk,'individual_prof_convolved_kcarta_AIRS_crisHI_*.mat'));
   flist = dir(strcat(dpk,'convolved_kcarta*radiances.mat'));
   for i=1:704 names{i}=flist(i).name; end
   fnum = 704;
end
[sorted_names,ndx] = natsortfiles(names);
clear krc kra;
for ip = 1:fnum
  load(strcat(dpk,sorted_names{ip}));
  for jp=1:6 krc(:,jp,ip) = rcris_all(:,jp); end
  for jp=1:6 kra(:,jp,ip) = rairs_all(:,jp); end
  fprintf('.');
end
fprintf('\n');
% This one off using regr49_1013_420ppm_2235g4.op.rtp (this set only nadir) to be
%  compared w/ sar_oct16_r49_1013_420_6a.rtp
fpk = '/asl/s1/chepplew/projects/sarta/cris_hr/xconvolved_kcarta_AIRS_crisHI.mat';
KJ  = load(fpk);
krc = KJ.rcris_all;

if(strcmp(csens,'CRIS'))
  freq  = fcris;
  bkc   = real(rad2bt(fcris,krc(:,iang,1:fnum-1)));
  bkcm  = nanmean(bkc,3);
  bsc   = real(rad2bt(hds.vchan,pds.rcalc(:,1:fnum-1)));
  bscm  = nanmean(bsc,2);
  bbias = bscm - bkcm(:,1);
  bias51 = bsc(51,:) - squeeze(bkc(51,1,:))';
  bias102 = bsc(102,:) - squeeze(bkc(102,1,:))';     (102:710.6250 cm-1)
  bias406 = bsc(406,:) - bkc(406,:);
  btstd = nanstd(squeeze(bkc(:,1,:)) - bsc,0,2);
%  bsstd = nanstd(bsc,1,2);
%  junk  = single(hamm_app(double(pds.robs1)));
%  bso   = real(rad2bt(hds.vchan,junk));
%  bsom  = nanmean(bso,2);
end
if(strcmp(csens,'AIRS'))
  bka   = real(rad2bt(fairs,kra(:,:,1:fnum-1)));
  bkam  = nanmean(bka,3);
  bbias = bscm - bkam(:,1);
  btstd = nanstd(squeeze(bka(:,iang,:)) - bsc,0,2);
  freq  = fairs;
end

%{
 figure(1);clf;h1=subplot(2,1,1);plot(fcris,bkcm(:,1),'-',fcris,bscm,'-');axis([640 1100 215 290]);
  title('SAF704 kcarta vs sarta mean (u), bias (l)');ylabel('BT K');
  h2=subplot(2,1,2);plot(fcris,bbias,'-',fcris,btstd,'-');axis([640 1100 -0.5 1]);grid on;
  legend('bias','std.dev','Location','best');xlabel('wavenumber cm^{-1}');

addpath /asl/matlib/rtptools
  mm_water = mmwater_rtp(hdr,pdr);
 figure(3);clf;scatter(mm_water(1:end-1),bias406,[],abs(pdr.rlat(1:end-1)));colorbar     
 figure(3);clf;scatter(abs(pdr.rlat(1:end-1)),bias406,[],mm_water(1:end-1));colorbar  
%}
%   ----------------- %%%%%%%%%% -----------------------------------------
%  2:             SAF704 sets (use zeroemiss)
%                     %%%%%%%%%%
dpk = ['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES'...
       '/JUNK_SAF_704_profiles_25-May-2016_xmb_zeroemiss/']; 
kcname = 'convolved_kcarta_crisHI_LBLRTM.mat';
load(strcat(dpk,kcname));
bkc   = real(rad2bt(fcris,rcris_all));
bkcm  = nanmean(bkc,2);
%
dpk = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/RAD1100/';
kcname = 'convolved_kcarta_RAD1100_*_radiances.mat';
flist = dir(strcat(dpk,kcname));
ropro = [100:109 10 110:119 11 120:129 12 130:139 13 140:149 14 ...
         150:159 15 160:169 16 170:179 17 180:189 18 190:199 19 1 ...
	 200:209 20 210:219 21 220:229 22 230:239 23 240:249 24 ...
	 250:259 25 260:269 26 270:279 27 280:289 28 290:299 29 2 ...
	 300:309 30 310:319 31 320:329 32 330:339 33 340:349 34 ...
	 350:359 35 360:369 36 370:379 37 380:389 38 390:399 39 3 ...
	 400:409 40 410:419 41 420:429 42 430:439 43 440:449 44 ...
	 450:459 45 460:469 46 470:479 47 480:489 48 490:499 49 4 ...
	 500:509 50 510:519 51 520:529 52 530:539 53 540:549 54 ...
	 550:559 55 560:569 56 570:579 57 580:589 58 590:599 59 5 ...
	 600:609 60 610:619 61 620:629 62 630:639 63 640:649 64 ...
	 650:659 65 660:669 66 670:679 67 680:689 68 690:699 69 6 ...
	 700:704 70:79 7 80:89 8 90:99 9];   % the order the mat files are loaded.
clear krc kra;
for ip = 1:numel(ropro)
  load(strcat(dpk,flist(ip).name));
  for jp=1:6 krc(:,jp,ropro(ip)) = rcris_all(:,jp); end
  for jp=1:6 kra(:,jp,ropro(ip)) = rairs_all(:,jp); end
end
%
bkc   = real(rad2bt(fcris,krc(:,:,1:704)));
bkcm  = nanmean(bkc,3);
%
sarfile = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_SAF704_2235.rtp';
sarfile = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_REGR_SAF704_25May2016_xmb_zeromiss.rtp';
[fpath fname fext] = fileparts(sarfile);
[hds has pds pas]  = rtpread(sarfile);
bsc   = real(rad2bt(hds.vchan,pds.rcalc(:,:)));
bscm  = nanmean(bsc,2);

whos bsc bscm bkc bkcm bbias btstd
%
bbias = bscm - bkcm(:,1);
btstd = nanstd(squeeze(bkc(:,:)) - bsc,0,2);

% -------------------------- END of SAF704 ----------------------------------

whos bsc bkc bka btstd bscm bkcm bkam bbias;

% plotting section
txtbands = {'LW','MW','SW'};
cbands   = [640, 1100; 1200, 1750; 2150, 2560];
abands   = [640, 1150; 1210, 1650; 2150, 2650];
%figure(1);clf;plot(fcris,bkcm(:,1),'.-',hds.vchan,bscm,'.-');grid on; 
%  axis([640 1100 215 280]);
if(strcmp(sens,'CRIS')) 
   bands = cbands; bkxm = bkcm; 
   pname = ['cris_r49_sarta_full_th2_vs_kcarta_ang' sprintf('%02.0f',scang(iang))];
   ctitle = ['sarta cris_oct16 vs kcarta r49 pCO2x1.05 ang ' sprintf('%02.0f',scang(iang))];
end
if(strcmp(sens,'AIRS')) 
   bands = abands; bkxm = bkam; 
   pname = 'airs_r49_sarta_vs_kcarta';
   ctitle = 'AIRS sarta vs kcarta regr 49';
end
for ib=1:3
figure(ib);clf;h1=subplot(2,1,1);plot(freq,bkxm(:,1),'.-',hds.vchan,bscm,'.-');grid on;
  ylabel('BT (K)');
  axis([bands(ib,:) 215 285]);ht=title(ctitle);set(ht,'interpreter','none');
  h2=subplot(2,1,2);plot(freq,bbias, '-',freq,btstd,'-');xlabel('wavenumber cm-1');
  grid on;axis([bands(ib,:) -0.5 1.0]);legend('bias','stdv','Location','northEast');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]) 
   %aslprint(['./figs/' pname '_' txtbands{ib} '.png']);
   %saveas(gcf,['./figs/' pname '_' txtbands{ib} '.png']);
end
%{
  saveas(gcf,'./figs/sarta_kc_regr49_LW_optr_v5.png','png');
%}
% find(fcris>1240,1) = 775 fcris(1832)=2300 wn
figure(2);clf;h1=subplot(2,1,1);plot(fcris,btstd,'-',fcris,bkcm - bscm,'-');
  grid on;axis([bands(2,:) -0.1 0.4]);

h2=subplot(2,1,2);plot([1:48],bkc(1832,:) - bsc(1832,:),'o-');grid on;

% -----------------------------------------------
% 3: Compare CALC with OBS from CCAST granule
% -----------------------------------------------
addpath /asl/packages/airs_decon/source      % hamm_app.m
ocfile = '/asl/s1/chepplew/projects/sarta/cris_hr/cris_2015d151_sar_wcon8_gran.rtp';

[hdc hac pdc pac] = rtpread(ocfile);
boc   = real(rad2bt(hdc.vchan,pdc.rcalc));
bocm  = nanmean(boc,2);
junk  = single(hamm_app(double(pdc.robs1)));
boo   = real(rad2bt(hdc.vchan,junk));
boom  = nanmean(boo,2);

ocbias = bocm - boom;
ocstd  =  nanstd(boc - boo,1,2);

figure(2);clf;h1=subplot(2,1,1);plot(hdc.vchan,bocm,'-',hdc.vchan,boom,'-');
  grid on;axis([bands(1,:) 210 300]);
  h2=subplot(2,1,2);plot(hdc.vchan,ocbias,'-',hdc.vchan,ocstd,'-');
  grid on;axis([bands(1,:) -1.2 1.2]);
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]) 

% ----------------------------------------------
% Compare merged coefficient sets for CrIS hrg4.
% ----------------------------------------------
addpath /home/chepplew/projects/sarta/matlabcode

[ichraw1 fchraw1 coefraw1 inforaw1] = rdcoef(1,0,'/home/chepplew/gitLib/ftc_dev/run/merg_set1_fow_coef.dat');
[ichraw2 fchraw2 coefraw2 inforaw2] = rdcoef(2,0,'/home/chepplew/gitLib/ftc_dev/run/merg_set2_fwo_coef.dat');
figure(3);clf;plot(ichraw1,coefraw1(:,40,16),'.-',ichraw2,coefraw2(:,40,26),'o');grid on;xlim([200 230]);
cut1=textread('/home/chepplew/gitLib/ftc_dev/chanLists/list_crisSet1_hrg4');
cut2=textread('/home/chepplew/gitLib/ftc_dev/chanLists/list_crisSet2_hrg4');
figure(3);clf;plot(cut1(:,1),cut1(:,2),'o',cut2(:,1),cut2(:,2),'o');grid on;

FNC1 = '/home/chepplew/gitLib/ftc_dev/run/cut_set1_fow_coef.dat';
FNC2 = '/home/chepplew/gitLib/ftc_dev/run/cut_set2_fwo_coef.dat';
FNC1 = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/set1_hrg4.dat';
FNC2 = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/set2_hrg4.dat';
FNC3 = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/set3_hrg4.dat';
FNCO = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/optran_hrg2.dat';
[ichcut1 fchcut1 coefcut1 infocut1] = rdcoef(1,0,FNC1);
[ichcut2 fchcut2 coefcut2 infocut2] = rdcoef(2,0,FNC2);
[ichcut3 fchcut3 coefcut3 infocut3] = rdcoef(3,0,FNC3);
[ichopt  fchopt  coefopt  infoopt]  = rdcoef(9,1,FNCO);
figure(2);clf;plot(ichcut1,fchcut1,'.',ichcut2,fchcut2,'o');grid on;
figure(3);clf;plot(ichcut1,coefcut1(:,40,16),'.-',ichcut2,coefcut2(:,40,26),'o');grid on; 
  xlim([180 230]);

% 20 June 2016 : compare methane coefficients ~ 1300 wn. (IASI vs mine)
FNC3S = '/asl/data/sarta_database/Data_IASI_may09/Coef/set3.dat';
FNC3C = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/set3_hrg4.dat';
[ich3s fch3s coef3s info3s] = rdcoef(3,0,FNC3S);
[ich3c fch3c coef3c inof3c] = rdcoef(3,0,FNC3C);

whos ich3s fch3s coef3s ich3c fch3c coef3c
figure(1);clf;plot(fch3s,coef3s(:,40,1),'-',fch3c,coef3c(:,40,1),'-');grid on;xlim([1285 1325]);

% 20 June 2016 : compare optran coefficients ~ 1300 wn (IASI vs mine)
FNCOptS = '/asl/data/sarta_database/Data_IASI_may09/Coef/optran.dat';
FNCOptC = '/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_may16/Coef/optran_hrg4.dat';
[ichOps fchOps coefOps infoOps] = rdcoef(9,1,FNCOptS);
[ichOpc fchOpc coefOpc infoOpc] = rdcoef(9,1,FNCOptC);
whos ichOps fchOps coefOps ichOpc fchOpc coefOpc
figure(3);clf;plot(fchOps,coefOps(:,40,1),'-',fchOpc,coefOpc(:,40,1),'-');grid on;xlim([1285 1325]);

%%%%%%%%%%%%%%%%%%%%%%%% Test with two surfaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
% 4:   49 Regression profiles at 6 angles and 2 surfaces
% ----------------------------------------------------------------------
% Data are arranged in order from  angles 0, 28.6101, 38.5184, 45.2214, 46.3539, 48.3250
% for 49 profiles each:  294 profiles, (6 angs, sea)
% original: ['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/'...
%       'regr_rtp_6angs_49profs_2surfaces.rtp'];
fnrtp = [/asl/s1/chepplew/projects/sarta/cris_hr/'...
         'regr_rtp_6angs_49profs_1013mb_seaemis_2235.rtp'];  
         'regr_rtp_6angs_49profs_1013mb_unitemis_2235.rtp'];
[hdr har pdr par] = rtpread(fnrtp);

XJ=load(['/asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/'...
         'JUNK/xconvolved_kcarta_crisHI_UMBCODs_g2_3_6_LBLRTMODs.mat']);
	 'JUNK/xconvolved_kcarta_crisHI_xconvolved_kcarta_crisHI_1013mb_unitemiss.mat']);
	 'JUNK/xconvolved_kcarta_crisHI_xconvolved_kcarta_crisHI_1013mb_seaemiss.mat']);
% fcris [2235x1] rcris_all [2235x588] satzen,scanang,solzen,stemp [1x588], rtp_file:
% select view angle
satzen = unique(XJ.satzen);
scang  = unique(XJ.scanang);
iang   = 1; sfc = 1;
indx   = [iang*49 - 48:iang*49];
freq   = XJ.fcris;
bkc    = real(rad2bt(XJ.fcris,XJ.rcris_all(:,indx)));    % 
bkcm   = nanmean(bkc,2);

sarfile = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_test2.rtp';
sarfile = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_kc_r49_1013_seaemiss.rtp';
[hds has pds pas] = rtpread(sarfile);
bsc   = real(rad2bt(hds.vchan, pds.rcalc(:,indx)));      % 
bscm  = nanmean(bsc,2);
whos bkc bkcm bsc bscm
bbias = bscm - bkcm;
btstd = nanstd(bkc - bsc,0,2);


% -----------------------------------------------------------------------
% convert pressure on levels to pressure of layer between the two levels.

pN = pdr.plevs(1:100,:)-pdr.plevs(2:101,:);
pD = log(pdr.plevs(1:100,:) ./ pdr.plevs(2:101,:));
pdr.plays = zeros(size(pdr.plevs));
pdr.plays(1:100,:) = pN ./ pD;

% -----------------------------------------------------------------------
