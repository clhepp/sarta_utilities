% perturb49_kc_vs_sarta.m

% comparison of TOA radiances derived from kcarta and SARTA for hi-res CrIS
% with minor gas profiles perturbed as follows:
% profiles: 
  pind =  {{'unp',  0,  10,20,30,40,50,60},...
           {'WV',   1,  01,11,21,31,41,51},... 
           {'CO2',  2,  02,12,22,32,42,52},...
           {'O3',   3,  03,13,23,33,43,53},...
           {'N2O',  4,  04,14,24,34,44,54},...
           {'CO',   5,  05,15,25,35,45,55},...
           {'CH4',  6,  06,16,26,36,46,56},...
           {'SO2',  9,  07,17,27,37,47,57},... 
           {'HNO3',12,  08,18,28,38,48,58},...
           {'T',   99,  09,19,29,39,49,59}};
%
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
warning 'off'
%
% Perturbed profiles: 
%original: /asl/s1/sergio/home/MATLABCODE/REGR_PROFILES/RUN_KCARTA/testperturb.rp.rtp
fnrtp = '/home/chepplew/data/sarta/prod_2016/sarta_data/cris_hr/testperturb_2235.rp.rtp';
fnrtp = '/home/chepplew/data/sarta/prod_2018/sarta_data/r49_perturb_1013m_400p_1e_2235g4.rtp';

% kcarta TOA radiances:
% original: /asl/s1/sergio/home/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/convolved_kcarta_crisHI.mat
KP = load('/asl/s1/chepplew/projects/sarta/cris_hr/convolved_kcarta_crisHI.mat');
KP = load('/home/chepplew/data/sarta/prod_2016/sarta_data/cris_hr/convolved_kcarta_crisHI.mat');

kpath='/home/chepplew/projects/kcarta/run/JUNK/';
fnlst = dir([kpath 'individual_prof_convolved_kcarta_AIRS_crisHI_*.mat']);
% reorder these in profile number order
for i=1:60 
  fnparts{i} = strsplit(fnlst(i).name, {'_','.'}); 
  prfnums(i) = str2num(cell2mat(fnparts{i}(7)));
end
  [B IB] = sort(prfnums);  
  krc = []; 
for i=IB 
  x   = load([kpath fnlst(i).name]);
  krc = [krc x.rcris_all];                 % [2235 x 60] 
end

% SARTA TOA radiances:
spath = '/home/chepplew/data/sarta/prod_2018/sarta_data/';
fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_400_7set_testperturb.rtp';
fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_400p_7s_so2_hno3_n2o_tp.rtp';
fnsar = '/home/chepplew/data/sarta/prod_2016/sarta_data/cris_hr/sar_400p_7s_co2_opt_nte_tra_tp.rtp';
fnsar = [spath 'sar_crisg4_mar18_basic_r49_tpert.rtp'];

[hds has pds pas] = rtpread(fnsar);

% align the two spectral grids
[xf xi xj] = intersect(x.fcris, hds.vchan);

% get the un-perturbed profiles, from sarta (s) and kcarta (k)
inUN  = cell2mat(pind{1}(3:8));
bsun  = real(rad2bt(hds.vchan, pds.rcalc(:,inUN)));
bsunm = nanmean(bsun,2);
%bkun  = real(rad2bt(KP.fcris, KP.rcris_all(:,inUN)));
bkun  = real(rad2bt(x.fcris, krc(:,inUN)) );
bkunm = nanmean(bkun,2);

bbias = bkunm(xi) - bsunm(xj);
figure(2);clf;plot(hds.vchan(1:2211), bbias(1:2211),'-');

% get the CO2 perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing  ' pind{3}(1)]);  ip=3;
inCO2  = cell2mat(pind{3}(3:8));
bsco2  = real(rad2bt(hds.vchan, pds.rcalc(:,inCO2)));
bsco2m = nanmean(bsco2,2);
bkco2  = real(rad2bt(KP.fcris, KP.rcris_all(:,inCO2)));
bkco2  = real(rad2bt(x.fcris, krc(:,inCO2)) );
bkco2m = nanmean(bkco2,2);

bkbias = bkunm-bkco2m;
bkstd  = nanstd(bkun - bkco2,1,2);
bsbias = bsunm-bsco2m;
bsstd  = nanstd(bsun - bsco2,1,2);

% get the O3 perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing  ' pind{4}(1)]);  ip=4;
inO3  = cell2mat(pind{4}(3:8));
bso3  = real(rad2bt(hds.vchan, pds.rcalc(:,inO3)));
bso3m = nanmean(bso3,2);
bko3  = real(rad2bt(KP.fcris, KP.rcris_all(:,inO3)));
bko3  = real(rad2bt(x.fcris, krc(:,inO3)) );
bko3m = nanmean(bko3,2);

bkbias = bkunm-bko3m;
bkstd  = nanstd(bkun - bko3,1,2);
bsbias = bsunm-bso3m;
bsstd  = nanstd(bsun - bso3,1,2);

% get the N2O perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing ' pind{5}(1)]);  ip=5;
inN2O  = cell2mat(pind{5}(3:8));
bsn2o  = real(rad2bt(hds.vchan, pds.rcalc(:,inN2O)));
bsn2om = nanmean(bsn2o,2);
bkn2o  = real(rad2bt(KP.fcris, KP.rcris_all(:,inN2O)));
bkn2o  = real(rad2bt(x.fcris, krc(:,inN2O)) );
bkn2om = nanmean(bkn2o,2);

bkbias = bkunm-bkn2om;
bkstd  = nanstd(bkun - bkn2o,1,2);
bsbias = bsunm-bsn2om;
bsstd  = nanstd(bsun - bsn2o,1,2);

% get the CO perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing ' pind{6}(1)]);  ip=6;
inCO  = cell2mat(pind{6}(3:8));
bsco  = real(rad2bt(hds.vchan, pds.rcalc(:,inCO)));
bscom = nanmean(bsco,2);
bkco  = real(rad2bt(KP.fcris, KP.rcris_all(:,inCO)));
bkco  = real(rad2bt(x.fcris, krc(:,inCO)) );
bkcom = nanmean(bkco,2);

bkbias = bkunm-bkcom;
bkstd  = nanstd(bkun - bkco,1,2);
bsbias = bsunm-bscom;
bsstd  = nanstd(bsun - bsco,1,2);

% get the CH4 perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing ' pind{7}(1)]);  ip=7;
inCH  = cell2mat(pind{7}(3:8));
bsch  = real(rad2bt(hds.vchan, pds.rcalc(:,inCH)));
bschm = nanmean(bsch,2);
bkch  = real(rad2bt(KP.fcris, KP.rcris_all(:,inCH)));
bkch  = real(rad2bt(x.fcris, krc(:,inCH)) );
bkchm = nanmean(bkch,2);

bkbias = bkunm-bkchm;
bkstd  = nanstd(bkun - bkch,1,2);
bsbias = bsunm-bschm;
bsstd  = nanstd(bsun - bsch,1,2);

% get the SO2 perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing ' pind{8}(1)]);  ip=8;
inSO  = cell2mat(pind{8}(3:8));
bsso  = real(rad2bt(hds.vchan, pds.rcalc(:,inSO)));
bssom = nanmean(bsso,2);
%bkso  = real(rad2bt(KP.fcris, KP.rcris_all(:,inSO)));
bkso  = real(rad2bt(x.fcris, krc(:,inSO)) );
bksom = nanmean(bkso,2);

bkbias = bkunm-bksom;
bkstd  = nanstd(bkun - bkso,1,2);
bsbias = bsunm-bssom;
bsstd  = nanstd(bsun - bsso,1,2);

bbias = bkunm(xi) - bsunm(xj);
figure(2);clf;plot(x.fcris,bkbias,'-', hds.vchan(1:2211),bsbias(1:2211),'-');
  grid on; xlim([1250 1450]);legend('kcarta','sarta');
  xlabel('wavenumber cm^{-1}');ylabel('nominal minus perturbed (K)')
  title('kcarta vs sarta calc. 10% SO2 perturb');
  
% get the HNO3 perturbed profiles, from sarta (s) and kcarta (k)
disp(['doing ' pind{9}(1)]);  ip=9;
inHN  = cell2mat(pind{9}(3:8));
bshn  = real(rad2bt(hds.vchan, pds.rcalc(:,inHN)));
bshnm = nanmean(bshn,2);
bkhn  = real(rad2bt(KP.fcris, KP.rcris_all(:,inHN)));
bkhnm = nanmean(bkhn,2);

bkbias = bkunm-bkhnm;
bkstd  = nanstd(bkun - bkhn,1,2);
bsbias = bsunm-bshnm;
bsstd  = nanstd(bsun - bshn,1,2);

%{
addpath /asl/matlib/plotutils
txtbands = {'LW','MW','SW'};
cbands   = [640, 1100; 1200, 1750; 2150, 2560; 2150, 2250];

for ib = 1:3
figure(ib+3);clf;h1=subplot(2,1,1);plot(KP.fcris, bkbias,'-', KP.fcris,bkstd,'--');grid on;
  axis([cbands(ib,:) -Inf Inf]); title(['kcarta v SARTA perturb ' char(pind{ip}(1)) ' bias, stddev']); 
  ylabel('BT (K)');legend('signal','std.dev','Location','North');
  h2=subplot(2,1,2);plot(hds.vchan, bsbias,'-', hds.vchan, bsstd,'--');grid on;
  axis([cbands(ib,:) -Inf Inf]);xlabel('wavenumber cm-1');ylabel('BT (K)');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1]) 
  %aslprint(['./figs/kc_vs_sar_full_p' pind{ip}{1} '_' txtbands{ib} '_290916a.png']);  
end  
  
% ---------------------------------------
% Section 2 - dedicated pCO2 files
% ---------------------------------------
dpk = {'/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/RAD1013_CO2x0.95/',...
       '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/RAD1013_CO2x1.0/',...
       '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/RAD1013_CO2x1.05/'};
ropro = [10:19 1 20:29 2 30:39 3 40:49 4:9];   % the order the mat files are loaded.
iang = 1;
for jc = 1:3
  clear krc kra;
  for ip = 1:49
    load(strcat(dpk{jc},flist(ip).name));
    for jp=1:6 krc(:,jp,ropro(ip)) = rcris_all(:,jp); end
    for jp=1:6 kra(:,jp,ropro(ip)) = rairs_all(:,jp); end
  end

  bkc{jc}   = real(rad2bt( fcris, squeeze(krc(:,iang,:)) ));     % nadir
  bkcnm{jc} = nanmean(bkc{jc},2);
end
% Alternate set of kcarta files
dpk = {'/asl/s1/chepplew/projects/sarta/cris_hr/xconvolved_kcarta_regr49_1013_400ppm_nadir_crisHI.mat',...
       '/asl/s1/chepplew/projects/sarta/cris_hr/xconvolved_kcarta_regr49_1013_410ppm_nadir_crisHI.mat',...
       '/asl/s1/chepplew/projects/sarta/cris_hr/xconvolved_kcarta_regr49_1013_420ppm_nadir_crisHI.mat'};
clear bkc bkcnm;
for jc=1:3
  KJ=load(dpk{jc});
  bkc{jc}   = real(rad2bt( KJ.fcris, KJ.rcris_all ));
  bkcnm{jc} = nanmean(bkc{jc},2);
end
%{
   figure(4);clf;hold on; for jc=1:3 plot(fcris,bkcnm{jc}); end;       
   figure(4);clf;plot(fcris,bkcnm{2} - bkcnm{1},'-', fcris, bkcnm{2} - bkcnm{3},'-')       
   xlim([640 1100]);grid on;legend('410-400ppm','410-400ppm','Location','north');       
   title('kcarta pCO2 relative to 400ppm');xlabel('wn cm^{-1}');ylabel('BT K');
   %saveas(gcf,'./figs/kcarta_pCO2_supp_LW.png','png')
%}
   
%fnsar = {'/asl/s1/chepplew/projects/sarta/cris_hr/sar_oct16_r49_1013_380_6a.rtp',
fnsar =  {'/asl/s1/chepplew/projects/sarta/cris_hr/sar_oct16_r49_1013_400_6a.rtp',...
	 '/asl/s1/chepplew/projects/sarta/cris_hr/sar_oct16_r49_1013_410_6a.rtp',...
	 '/asl/s1/chepplew/projects/sarta/cris_hr/sar_oct16_r49_1013_420_6a.rtp'};
	 
clear bsc bscnm;
for jc = 1:3;
  [hds has pds pas] = rtpread(fnsar{jc});
  bsc{jc}   = real(rad2bt(hds.vchan,pds.rcalc));     % nadir only
  bscnm{jc} = nanmean(bsc{jc},2);
end
fc = hds.vchan;

%{
   figure(5);clf;hold on; for jc=1:3 plot(fcris,bscnm{jc}); end;       
   figure(5);clf;plot(fc,bscnm{2} - bscnm{1},'-',fc,bscnm{2} - bscnm{3},'-');
   xlim([640 1100]);grid on;legend('410-400ppm','410-420ppm','Location','north');       
   title('sarta.oct16 pCO2 relative to 400ppm');xlabel('wn cm^{-1}');ylabel('BT K');
   %saveas(gcf,'./figs/sar_oct16_pCO2_supp_LW.png','png')

% direct comparison with kcarta
  dsar_410  = (bscnm{2} - bscnm{1});   dkc_410 = (bkcnm{2} - bkcnm{1});
  dbias_410 =  dsar_410 - dkc_410;
  figure(6);clf;plot(fc,dsar_410,'-',fc,dkc_410,'-',fc,dbias_410,'-');
    xlim([640 1100]);grid on;title('sarta vs kcarta pCO2 deltas');xlabel('wn cm^{-1}');
    ylabel('BT K');legend('sar 410-400','kc 410-400','dsar-dkc');
   % saveas(gcf,'./figs/sar_kc_pCO2_410-400ppm_deltas_LW.png','png');
%}
 
% -------------------------------------------
% Setting up the RTP files for this pCO2 test
% -------------------------------------------
% get hdr2 with channel vctors for CrIS hiRes w/ 4 guard channels
load('rtp_head_str_2235g4.mat')

%
fnrtp='/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1013_400ppm.op.rtp';
[hdr har pdr par] = rtpread(fnrtp);
pdr_orig = pdr;

co2_0p95 = pds.gas_2*0.95;
pdr2 = pdr_orig;
pdr2.gas_2 = co2_0p95;
fortp = 'regr49_1013_380ppm_2235g4.op.rtp';
rtpwrite(fortp,hdr2,har,pdr2,par);

co2_1p05 = pdr_orig.gas_2*1.05;
pdr2 = pdr_orig;
pdr2.gas_2 = co2_1p05;                   
fortp = 'regr49_1013_420ppm_2235g4.op.rtp';
rtpwrite(fortp,hdr2,har,pdr2,par); 
fortp = 'regr49_1013_400ppm_2235g4.op.rtp';
rtpwrite(fortp,hdr2,har,pdr_orig,par);

fnrtp='/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1013_400ppm.op.rtp';
[hdr har pdr par] = rtpread(fnrtp);
pdr_orig = pdr;
co2_410 = pdr.gas_2*1.025;
fortp = 'regr49_1013_410ppm.op.rtp';
pdr2=pdr_orig;
pdr2.gas_2 = co2_410;
rtpwrite(fortp,hdr,har,pdr2,par);
