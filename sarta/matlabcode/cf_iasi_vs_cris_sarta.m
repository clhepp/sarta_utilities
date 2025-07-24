% cf_iasi_vs_cris_sarta.m

% compare results of running IASI_SARTA deconvolved to CrIS, with CrIS_SARTA

% -------------------------------------------------------------
% Part 1:  using test atmospheric profiles (already on klayers)
% -------------------------------------------------------------
cd /home/chepplew/projects/sarta/cris_hr

addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil            % mktemp.m

bands = {[640 1100], [1200 1800], [2150 2560]};
plevs = textread('/home/sergio/MATLABCODE/airslevels.dat');

iasi_sarta = '/asl/packages/sartaV108/BinV201/sarta_iasi_may09_wcon_nte';

rtpfile = '/home/sergio/MATLABCODE/REGR_PROFILES/ECMWF_SAF_137Profiles/save_SAF_704_profiles_29-Apr-2016.op.rtp';

[head hattr prof pattr] = rtpread(rtpfile);

% need to populate head.vchan and ichan w/ IASI grid.
addpath /asl/packages/iasi_decon
iasi = iasi_params;
head.ichan = [1:8461]';
head.vchan = iasi.freq;
% for testing nonLTE need to hack the solzen value (all SAF profs have 135 deg)
prof.solzen = single(20.0 * ones(1,704));

 tmp = mktemp();
  outfiles = rtpwrite_12(tmp,head,hattr,prof,pattr);
  s1Path = '/asl/s1/chepplew/tmp/';

  ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
  ofn_1 = [tmp '.kla_1'];  ofn_2 = [tmp '.kla_2'];
  ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];


 % run sarta on first half
  eval(['! ' iasi_sarta ' fin=' ifn_1 ' fout=' ofn_3 ' > sartastdout1.txt']);

 % run sarta on second half
  eval(['! ' iasi_sarta ' fin=' ifn_2 ' fout=' ofn_4 ' > sartastdout1.txt']);

 % read the results files back in
  cfin = [tmp '.sar'];

  [hd ha pd pa] = rtpread_12(cfin);

unlink tmp; unlink ifn_1; unlink ifn_2; unlink ofn_3; unlink ofn_4;

% deconvolve rcalc to CrIS
addpath /asl/packages/iasi_decon 
addpath /asl/packages/ccast/source           % inst_params.m
addpath /asl/packages/airs_decon/source      % hamm_app.m

opt1 = struct;
opt1.hapod   = 1;
opt1.resmode = 'hires2';
opt1.nguard  = 2;

[crad cfreq] = iasi2cris(pd.rcalc, hd.vchan, opt1);

i2cbt  = real(rad2bt(cfreq, crad));
i2cbtm = nanmean(i2cbt,2);
  junk = single(hamm_app(double(pd.rcalc)));
ibt    = real(rad2bt(hd.vchan,junk));

figure(1);clf;plot(cfreq,i2cbt(:,1),'-',hd.vchan,ibt(:,1),'-');grid on;xlim(bands{3});

% -------------------------------------
% Now run CrIS_SARTA
% -------------------------------------
cris_sarta = '/home/chepplew/gitLib/sarta/sarta_cris_hrg2_mar16_wcon';
cris_sarta = '/home/chepplew/gitLib/sarta/sarta_cris_hrg2_mar16_wcon_nte';

% reload original rtpfile and decon to get hi-res CrIS grid
 rtpfile = '/home/sergio/MATLABCODE/REGR_PROFILES/ECMWF_SAF_137Profiles/save_SAF_704_profiles_29-Apr-2016.op.rtp';
[head hattr prof pattr] = rtpread(rtpfile);

head.vchan = cfreq;       % [2223x1] carried forward for hires2
head.ichan = [1:2223]';
head.nchan = 2223;
head.vcmin = min(cfreq);
head.vcmax = max(cfreq);
% for testing nonLTE need to hack the solzen value (all SAF profs have 135 deg)
prof.solzen = single(20.0 * ones(1,704));

tmp = mktemp();
rtpwrite(tmp,head,hattr,prof,pattr);

ofn_5 = './cris_hrg2_nte_sar.rtp';

eval(['! ' cris_sarta ' fin=' tmp ' fout=' ofn_5 ' > sartastdout2.txt']);

[crHd crHa crPd crPa] = rtpread(ofn_5);

unlink tmp;

crBt = real(rad2bt(cfreq, crPd.rcalc));
crBtm = nanmean(crBt,2);

bias  = i2cbtm - crBtm;
stdev = nanstd(i2cbt - crBt,1,2);

figure(2);clf;h1=subplot(2,1,1);plot(cfreq,i2cbtm - crBtm,'-');grid on;
  xlim(bands{3});ylim([-5 10]);title('iasi.sarta.decon.hires.cris vs hires.cris.sarta.nte');
  ylabel('BT K');legend('mean bias');
h2=subplot(2,1,2);plot(cfreq,stdev,'-');grid on; xlim(bands{3});xlabel('wn');ylabel('BT K');
  legend('std dev')
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
  %saveas(gcf,'./figs/cf_iasi_vs_cris_nte_sarta_704profs_SW.png','png'); 


[crHd1 crHa1 crPd1 crPa1] = rtpread('./cris_hrg2_sar.rtp');
crBt1  = real(rad2bt(cfreq, crPd1.rcalc));
crBtm1 = nanmean(crBt1,2);
bias1  = i2cbtm - crBtm1; 
stdev1 = nanstd(i2cbt - crBt1,1,2);
figure(3);clf;h1=subplot(2,1,1);plot(cfreq,i2cbtm - crBtm1,'-');grid on;
  xlim(bands{3});ylim([-2 2]);

% -------------------------------------------------------------
% Part 2:  using day of interest DOI 2016.01.20 (already on klayers)
% -------------------------------------------------------------
addpath /asl/packages/airs_decon/source        % hamm_app.m

fnrtp = '/asl/rtp/rtp_cris_ccast_hires/clear/2016/cris_hr_ecmwf_d20160120_clear_isarta.rtp';
[head hatt prof patt] = rtpread(fnrtp);
fcris   = head.vchan;
inn     = find(prof.solzen > 90);
junk    = single(hamm_app(double(prof.robs1(:,inn))));
bto_a   = real(rad2bt(head.vchan,junk));
btom_a  = nanmean(bto_a,2);
junk    = single(hamm_app(double(prof.rcalc(:,inn))));
btc_a   = real(rad2bt(head.vchan,junk));
btcm_a  = nanmean(btc_a,2);
bbias_a = nanmean(bto_a - btc_a,2);
bstdv_a = nanstd(bto_a - btc_a, 0,2);

fnrtp = '/asl/rtp/rtp_cris_ccast_hires/clear/2016/cris_hr_ecmwf_d20160120_clear_csarta.rtp';
[head hatt prof patt] = rtpread(fnrtp);
% subset for night (eclipsed)
inn     = find(prof.solzen > 90);
junk    = single(hamm_app(double(prof.robs1(:,inn))));
bto_b   = real(rad2bt(head.vchan,junk));
btom_b  = nanmean(bto_b,2);
btc_b   = real(rad2bt(head.vchan,prof.rcalc(:,inn)));
btcm_b  = nanmean(btc_b,2);
bbias_b = nanmean(bto_b - btc_b,2);
bstdv_b = nanstd(bto_b - btc_b,0,2);


%{
bands = [640 1100; 1200 1760; 2150 2550];
figure(1);clf;h1=subplot(2,1,1);plot(fcris,btom_a(:,1),'-',fcris,btcm_a,'-');xlim([2150 2550]);
  title('2012.01.20 clear subset Obs Calc iasi.sarta');ylabel('BT K');
  h2=subplot(2,1,2);plot(fcris,bbias_a,'-',fcris,bstdv_a,'-');axis([2150 2550 -3 6]);grid on;
  xlabel('wavenumber cm-1');ylabel('d(BT) K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
 % saveas(gcf,'./figs/crisHR_bias_SW_night_20120120_iasiSarta.png','png');
figure(2);clf;h1=subplot(2,1,1);plot(fcris,btom_b,'-',fcris,btcm_b,'-');xlim([2150 2550]);
  title('2012.01.20 clear subset Obs Calc cris.sarta');ylabel('BT K');
  h2=subplot(2,1,2);plot(fcris,bbias_b,'-',fcris,bstdv_b,'-');axis([2150 2550 -3 6]);grid on;
  xlabel('wavenumber cm-1');ylabel('d(BT) K');
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])  
 % saveas(gcf,'./figs/crisHR_bias_SW_night_20120120_crisSarta.png','png');
%}

% Now re-run the sarta for cris hi-Res.
sarta_exe = '';
