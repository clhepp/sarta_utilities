%script: compare_sarta_and rtp_versions()

%
% Compare different versions of SARTA using R49 or RTP subsets
% Context: Updating SARTA for AIRS L1C, IASI, CrIS RTP production.
%   Sample RTP sett: 2019d018
%
% C. Hepplewhite: May 2020
% v2 CLH Nov2022. prod_2022,jul22, HITRAN2020 builds
%

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools              % rtpwrite_12
addpath /asl/matlib/aslutil

% Choose sensor
csens = 'CRIS_HR';
all_sens = {'AIRS_L1C','CRIS_LR','CRIS_HR','IASI','CHIRP'};

% --------------------------------------------------------
% The sarta executables used then the rtp profiles to use:
% --------------------------------------------------------
% AIRS
airs.sar.clr1  = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
airs.sar.sct1  = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
%sarta.new.clr  = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_SItest';
airs.sar.clr2  = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_jpl_tunmlt';
airs.sar.sct2  = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod';
%
airs.rtp.clr1  = '/asl/rtp/rtp_airicrad_v6/clear/2019/era_airicrad_day018_clear.rtp';
airs.rtp.rnd1  = '/asl/rtp/rtp_airicrad_v6/random/2019/era_airicrad_day018_random.rtp';
airs.rtp.afv1  = '/home/sbuczko1/Work/chepplew/2019_oldrta/018/';
airs.rtp.clr2  = '/home/sbuczko1/Work/chepplew/rtp_clear_2019_018_may19_prod.rtp';
airs.rtp.rnd2  = '/home/sbuczko1/Work/chepplew/rtp_random_2019_018_may19_prod.rtp';
airs.rtp.afv2  = '/home/sbuczko1/Work/chepplew/2019/018/';

% CriS FSR
cris.sar.clr1  = '/asl/packages/sartaV108/BinV201/sarta_crisg4_nov09_wcon_nte';
cris.sar.clr2  = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16';
cris.sar.clr3  = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_prod';
cris.sar.clr4  = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev';

% IASI

% CHIRP

% -------------------------------
% Standard 49 regression profiles
% -------------------------------
r49.rtp    = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_unitemis_seaemis_7angs_night.rtp';
[hd0, ha0, pd0, pa0] = rtpread(r49.rtp);
idx = struct;
uemis = unique(pd0.emis(9,:));
for i = 1:length(uemis)
  idx.emis{i} = find(pd0.emis(9,:) == uemis(i));
end

% 
Switch csens
  case 'CRIS_HR'
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    hd0.vchan = freq;
    hd0.ichan = idchan;
    hd0.nchan = length(idchan);
    %hd0.pmax  = 1013;
    tmp = mktemp();
    rtpwrite(tmp,hd0,ha0,pd0,pa0);
    ifn = tmp;      ofn = [tmp '.sar'];
    SARTA_CMD = [cris.sar.clr3 ' fin=' ifn ' fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log'];
    tic
      system(SARTA_CMD);
    toc
    [hds,~,ptemp,~] = rtpread(ofn);
    cris.bsc.clr3   = rad2bt(hd0.vchan, ptemp.rcalc);

    SARTA_CMD = [cris.sar.clr4 ' fin=' ifn ' fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log'];
    tic
      system(SARTA_CMD);
    toc
    [hds,~,ptemp,~] = rtpread(ofn);
    cris.bsc.clr4   = rad2bt(hd0.vchan, ptemp.rcalc);
  
  case 'AIRS_L1C'
  
  case 'IASI'
  
end


% -----------------------------------------------------------------
% Load existing RTP results and convert to BT
rtp.old.afov_dir = dir([rtp.old.allfov 'allfov_era_airicard_d2019018_*_prod.rtp']);
rtp.new.afov_dir = dir([rtp.new.allfov 'allfov_era_airicard_d2019018_*_prod.rtp']);

% load and convert to BT
[hd ha pd pa] = rtpread(rtp.old.clear);
fa = hd.vchan;
btc.old.clear  = rad2bt(fa, pd.rclr);
bto.old.clear  = rad2bt(fa, pd.robs1);
lat.old.clear  = pd.rlat;
lon.old.clear  = pd.rlon;
solz.old.clear = pd.solzen;

[hd ha pd pa]  = rtpread(rtp.new.clear);
btc.new.clear  = rad2bt(fa, pd.rclr);
bto.new.clear  = rad2bt(fa, pd.robs1);
lat.new.clear  = pd.rlat;
lon.new.clear  = pd.rlon;
solz.new.clear = pd.solzen;

[hd ha pd pa]  = rtpread(rtp.old.rand);
btc.old.randm  = rad2bt(fa, pd.rclr);
btd.old.randm  = rad2bt(fa, pd.rcld);
bto.old.randm  = rad2bt(fa, pd.robs1);
lat.old.randm  = pd.rlat;
lon.old.randm  = pd.rlon;
solz.old.randm = pd.solzen;

[hd ha pd pa] = rtpread(rtp.new.rand);
btc.new.randm = rad2bt(fa, pd.rclr);
btd.new.randm = rad2bt(fa, pd.rcld);
bto.new.randm = rad2bt(fa, pd.robs1);
lat.new.randm = pd.rlat;
lon.new.randm = pd.rlon;
solz.new.randm = pd.solzen;

igran=10;
[hd ha pd pa] = rtpread([rtp.old.allfov rtp.old.afov_dir(igran).name]);
btc.old.allfov = rad2bt(fa, pd.rclr);
btd.old.allfov = rad2bt(fa, pd.rcld);
bto.old.allfov = rad2bt(fa, pd.robs1);
lat.old.allfov = pd.rlat;
lon.old.allfov = pd.rlon;
solz.old.allfov = pd.solzen;

[hd ha pd pa] = rtpread([rtp.new.allfov rtp.new.afov_dir(igran).name]);
btc.new.allfov = rad2bt(fa, pd.rclr);
btd.new.allfov = rad2bt(fa, pd.rcld);
bto.new.allfov = rad2bt(fa, pd.robs1);
lat.new.allfov = pd.rlat;
lon.new.allfov = pd.rlon;
solz.new.allfov = pd.solzen;



% ----------------------------------------------------------------------
% PART 2 - double diff with kCARTA TOA rads using 49 regression profiles
% ----------------------------------------------------------------------
cd /home/chepplew/data/sarta/

srcdr  = '/home/chepplew/data/sarta/prod_2019/generic/';
srcrtp = [srcdr 'r49_1100_98lev_400p_unitemis_seaemis_7angs_night.rtp']; % or _1013_
kpath  = [srcdr '/kcarta/REGR49/400p_1100mb_98lev_7angs_unit_seaemis/'];
%
srcdr  = '/home/chepplew/data/sarta/prod_2022/generic/';
srcrtp = [srcdr 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp']; % or _1013_
kpath  = [srcdr '/kcarta/r49_1013_98lev_400p_unitemis_seaemis_7angs_night/'];

% Load the RTP data in preparation for processing.
[head hatt prof patt] = rtpread(srcrtp);
if(isfield(prof,'rcalc')) disp('Removing rcalc');
  prof = rmfield(prof,'rcalc');
end


fin = '';
switch csens
  case 'AIRS_L1C'
    hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
    freq = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    %freq    = load('/home/chepplew/myLib/data/airs_f_l1c_lls.txt');
    %idchan  = int32([1:2645])';
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);

    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);

    SARTAEXE = sarta.old.clear;
    ifn = tmp;      ofn = [tmp '.sar'];
    eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    [~,~,ptemp,~] = rtpread(ofn);
    calc.old.rad = ptemp.rcalc;
    calc.old.bsc = rad2bt(head.vchan, ptemp.rcalc);
%
    SARTAEXE = sarta.new.clear;
    eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    [~,~,ptemp,~] = rtpread(ofn);
    calc.new.rad = ptemp.rcalc;
    calc.new.bsc = rad2bt(head.vchan, ptemp.rcalc);

    fairs = head.vchan;
    [~, ib]  = sort(fairs);
  case 'CRIS_HR'



% Get the subsetting (to be used for plotting)
subs.ang = unique(prof.satzen);
idx = struct;
for i=1:7
  idx.a{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) == 1);
  idx.b{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) < 1);
end
subs.names = {'scan angles x 7','2 x surface emissivity: unit and sea'};

% use first 6 view angles for each surface subset:
iwnt = [];
for j=1:6
  iwnt = [iwnt idx.a{j}];
end
iwnt  = sort(iwnt);
fsort = fairs(ib);
bias = rad2bt(head.vchan(ib), nanmean(calc.old.rad(ib,iwnt),2)) ...
     - rad2bt(head.vchan(ib), nanmean(calc.new.rad(ib,iwnt),2));
figure(2);clf;plot(fsort, bias ,'-')
  title('old vs new sarta, mean regr49, night')
  ylabel('old minus new (K)'); xlim([600 2700]); xlabel('wvn cm-1')

addpath ~/gitLib/asl_sno/source/
  [~, iis] = sort(fairs);
   xx = a2c_good_chans(sort(fairs));

iil1c = setdiff([1:2834], xx.fill);
  hold on; plot(fsort(iil1c), bias(iil1c) ,'g.')
  
figure(2);clf;plot(fsort(iil1c), bias(iil1c) ,'-');xlim([600 2700]); grid on;
    title('old vs new sarta, mean regr49, night')
    ylabel('old minus new (K)'); xlabel('wvn cm-1')
% saveas(gcf,[phome 'old_vs_new_airs_l1c_sarta_regr49_mean_diff.fig'],'fig')

% --------------------------------------------------------------------
% PART 3 - take a day of random RTP, rerun calcs using old + new sarta
% ---------------------------------------------------------------------
[hd ha pd pa]  = rtpread(rtp.old.rand);
hd.pfields = 5;
pd = rmfield(pd,'rclr');
pd = rmfield(pd,'rcld');
FIN = '/scratch/tmp.rtp';
rtpwrite(FIN, hd, ha, pd, pa);
SARTAEXE=sarta.old.clear
FOUT = '/scratch/tmp.sar';
eval(['! ' SARTAEXE ' fin=' FIN '  fout=' FOUT ' > /home/chepplew/logs/sarta/sar_out.log']);
[~,~,ptemp,~] = rtpread(ofn);

% ---------------------------------------------------------
% plotting section
% ---------------------------------------------------------
% RTA report November 2022 (for CrIS and CHIRP
% add NaN in spectra to create break in the plot lines
freq = [hd0.vchan(1:721); NaN; hd0.vchan(722:1594); NaN; hd0.vchan(1595:end)];
bsc3 = [cris.bsc.clr3(1:721,:); NaN(1,686); cris.bsc.clr3(722:1594,:); ...
        NaN(1,686); cris.bsc.clr3(1595:end,:)];
bsc4 = [cris.bsc.clr4(1:721,:); NaN(1,686); cris.bsc.clr4(722:1594,:); ...
        NaN(1,686); cris.bsc.clr4(1595:end,:)];
fn=1;
fh1=figure(fn); clf;
set(gcf,'Resize','Off'); set(gcf,'Position',fh1.Position+[0 0 0 210]);
sb1=subplot(311);
  plot(freq, mean(bsc4,2),'-', freq,mean(bsc3,2),'-')
  grid on; xlim([600 2600]);ylabel('BT (K)');
  legend('mean BT','location','north')
  title('CrIS SARTA v2018 vs v2020')
  %%%annotation('textbox',[0.1 0.8 0.1 0.1],'String','a','FitBoxToText','on');
sb2=subplot(312);
  plot(freq, mean(bsc4,2)-mean(bsc3,2),'-')
  grid on; xlim([600 2600]);ylabel('delta-BT (K)');
  legend('mean bias','location','south')
sb3=subplot(313);
  plot(freq, std(bsc4,0,2),'-', freq,std(bsc3,0,2),'-')
  grid on; xlim([600 2600]);ylabel('BT (K)');xlabel('wavenumber (cm^{-1})');
  legend('std.dev','location','north');

linkaxes([sb1 sb2 sb3],'x');set(sb1,'xticklabel','');set(sb2,'xticklabel','');
pp1=get(sb1,'position');set(sb1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
pp2=get(sb2,'position');set(sb2,'position',[pp2(1) pp2(2)-pp2(4)*0.05 pp2(3) pp2(4)*1.05])
pp2=get(sb2,'position');set(sb2,'position',[pp2(1) pp2(2)+pp2(4)*0.05 pp2(3) pp2(4)*1.05])
pp3=get(sb3,'position');set(sb3,'position',[pp3(1) pp3(2)+pp3(4)*0.1 pp3(3) pp3(4)*1.1])




%{
% Clear
figure(1);clf;plot(fa, nanmean(btc.old.clear,2),'-'); hold on;
   plot(fa, nanmean(btc.new.clear,2),'-');
clf(1); plot(fa, nanmean(btc.old.clear,2) - nanmean(btc.new.clear,2),'-') 
  grid on; ylabel('old minus new (K)')
  title('2019d018 AIRS.L1C clear RTP subset. rclr')

iid.old = find(solz.old.clear < 90);
iid.new = find(solz.new.clear < 90);
iin.old = find(solz.old.clear > 90);
iin.new = find(solz.new.clear > 90);

clf(2);plot(fa, nanmean(bto.old.clear,2) - nanmean(bto.new.clear,2),'-') 
  grid on; ylabel('old minus new (K)')
  title('2019d018 AIRS.L1C clear RTP subset. robs')
  
ach = find(fa>899,1); 
pdf.new.clear = histcounts(btc.new.clear(ach,:), [260:1:320]);
pdf.old.clear = histcounts(btc.old.clear(ach,:), [260:1:320]);
clf(2);plot([260.5:1:319.5], pdf.new.clear,'.-', [260.5:1:319.5], pdf.old.clear,'.-')

% Random. clear & cloudy
figure(1);clf;plot(fa, nanmean(btc.old.randm,2),'-'); hold on;
   plot(fa, nanmean(btc.new.randm,2),'-');
clf(1); plot(fa, nanmean(btc.old.randm,2) - nanmean(btc.new.randm,2),'-') 
  xlim([600 2700]);grid on; ylabel('old minus new (K)')
  title('2019d018 AIRS.L1C random RTP subset. rclr')

% cloudy
figure(1);clf;plot(fa, nanmean(btd.old.randm,2),'-'); hold on;
   plot(fa, nanmean(btd.new.randm,2),'-');
clf(1); plot(fa, nanmean(btd.old.randm,2) - nanmean(btd.new.randm,2),'-') 
  xlim([600 2700]);grid on; ylabel('old minus new (K)')
  title('2019d018 AIRS.L1C random RTP subset. rcld')
% Random cloudy PDF
iid.old = find(solz.old.randm < 90);
iid.new = find(solz.new.randm < 90);
iin.old = find(solz.old.randm > 90);
iin.new = find(solz.new.randm > 90);

clf(1);plot(fa, nanmean(bto.old.randm,2) - nanmean(bto.new.randm,2),'-') 

ach = find(fa>899,1); 
pdf.new.randm = histcounts(btd.new.randm(ach,:), [220:1:320]);
pdf.old.randm = histcounts(btd.old.randm(ach,:), [220:1:320]);
clf(2);plot([220.5:1:319.5], pdf.new.randm,'.-', [220.5:1:319.5], pdf.old.randm,'.-')

figure; simplemap(lat.old.randm, lon.old.randm, btd.old.randm(ach,:))
%}
