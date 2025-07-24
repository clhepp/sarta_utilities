% compare_kc_sar_r49_pbl_airslay_generic
%
% same as comapre_kc_sar_r49_generic except for double differencing
%      airs.plevs vs pbl.plevs versions
%
% INPUT: opts structure with fields:
%        csens    {'AIRS_L1C,'AIRS_PBL',CRIS_HR','CRIS_PBL'...}
%        prod:    {'2025'}
%        build:   {'jan2025a'}
%        regset:  {'r49,'saf704'}
%  kCARTA results: TEST1: PBL.    TEST2: AIRSLAY.
%

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

% -----------------------------------------------------------------
% Check input parameters are valid
if( ~all(isfield(opts,{'csens','prod','build','regset'}))
  error('Please enter 4 option values'); return;
end

% Check which sensor to compare w/kcarta {'IASI','AIRS','CRIS'}
csens = upper(opts.csens);
allsens = {'AIRS_L1C','AIRS_PBL','CRIS_LR','CRIS_HR','CRIS_PBL',...
           'IASI','CHIRP','CHIRP_PBL'};
if(~ismember(csens,allsens)) error('Invalid Sensor'); return; end

% check production
prod = lower(opts.prod);
if(~ismember(prod,{'2025'}))
  error('invalid production date');
  return;
end
prod_run = ['prod_' prod];

% check build date
build = lower(opts.build);
if(~ismember(build,{'jan2025a'}))
  error('invalid build');
  return;
end

% -----------------------------------------------------------------------
% initialization
sartaexe.ap = '/home/chepplew/gitLib/sarta/bin/airs_pbl_dev';                  % airs.l1c.PBL
sartaexe.aa = '/home/chepplew/gitLib/sarta/bin/airs_l1c_p2025jul22_std_dev';   % airs.l1c.airslay
sartaexe.cp = '/home/chepplew/gitLib/sarta/bin/cris_hr_pbl_dev';               % cris.hr.PBL
sartaexe.ca = '/home/chepplew/gitLib/sarta/bin/cris_hr_p2025jul22_reg_dev';    % cris.hr.airslay
sartaexe.cx = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev';      % cris.hr.airslay
%
sartaexe.hp = '/home/chepplew/gitLib/sarta/bin/chirp_pbl_p2025jan25a_dev';
sartaexe.ha = '/home/chepplew/gitLib/sarta/bin/chirp_p2022jul22_dev';
% sergios Jacobian version:
sartaexe.serg = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020';

% Temprary scratch directory
tempPath = mktemp();

% kCARTA TOA rad calcs home
% original: /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/Test_ArbPLEVS_PBL/TEST[12]_Jan29_2025/
src.dir  = ['/home/chepplew/data/sarta/' prod_run '/generic/'];

switch prod
  case '2025'
    kc.home = '/home/chepplew/data/sarta/prod_2025/generic/kcarta/';
    kc.list1   = dir([kc.home 'TEST1/individual_prof_convolved*.mat']);
    kc.list2   = dir([kc.home 'TEST2/individual_prof_convolved*.mat']);
    src.rtp1   = [src.dir 'r49_1013_400p_8x2x2_2834_pbl.rtp']; 
    src.rtp2   = [src.dir 'r49_1013_400p_8x2x2_2834_airslay.rtp'];
end

% ---------------------------------------------------------------------
% set up rtp layers file for chosen sensor and command string
% ---------------------------------------------------------------------
switch csens
  case {'AIRS_L1C','AIRS_PBL'}
   ichan = hdfread('/home/chepplew/myLib/data/airs_l1c_srf_tables_lls_new.hdf','chanid');
   vchan = hdfread('/home/chepplew/myLib/data/airs_l1c_srf_tables_lls_new.hdf','freq');

  case {'CRIS_HR','CRIS_PBL'}
   load('/home/chepplew/myLib/data/fcris_hires_4grd.mat','ichan');;
   load('/home/chepplew/myLib/data/fcris_hires_4grd.mat','vchan');
  
  case {'CHIRP','CHIRP_PBL'}   % alt. chirp_1702_freq_4grd.mat
   %load('/home/chepplew/myLib/data/chirp_1702_freq.mat','ichan');;
   ichan = [1:1702];
   load('/home/chepplew/myLib/data/chirp_1702_freq.mat','vchan');
end
% load up the layers file and swap header spectrum
[head,hatt,prof,patt]  = rtpread(src.rtp1);    % 1: PBL
[head,hatt,prof,patt]  = rtpread(src.rtp2);    % 2: airslay
head.vchan = single(vchan(:));
head.ichan = int32(ichan(:));
head.nchan = length(ichan);
% write this LAYERS rtp file to scratch area ready for SARTA
fn.rtp_op1  = [tempPath '_' lower(csens) '_op1.rtp'];
   rtpwrite(fn.rtp_op1, head,hatt,prof,patt);
fn.rtp_sar1 = [tempPath '_' lower(csens) '_sar1.rtp'];
fn.rtp_op2  = [tempPath '_' lower(csens) '_op2.rtp'];
   rtpwrite(fn.rtp_op2, head,hatt,prof,patt);
fn.rtp_sar2 = [tempPath '_' lower(csens) '_sar2.rtp'];

switch csens
  case {'AIRS','AIRS_PBL'}
% AIRS : choose which SARTA exec to use (1: PBL, 2: airslay)
    command1 = [sartaexe.ap ' fin=' fn.rtp_op1 ' fout=' fn.rtp_sar1 ...
          ' > /home/chepplew/data/scratch/sar_out.log'];
    command2 = [sartaexe.aa ' fin=' fn.rtp_op2 ' fout=' fn.rtp_sar2 ...
          ' > /home/chepplew/data/scratch/sar_out.log'];
  case {'CRIS_HR','CRIS_HR_PBL'}
% CRIS
    command1 = [sartaexe.cp ' fin=' fn.rtp_op1 ' fout=' fn.rtp_sar1 ...
           ' > /home/chepplew/data/scratch/sar_out.log'];
    command2 = [sartaexe.ca ' fin=' fn.rtp_op2 ' fout=' fn.rtp_sar2 ...
          ' > /home/chepplew/data/scratch/sar_out.log'];
  case {'CHIRP','CIRP_PBL'}
% CHIRP
    command1 = [sartaexe.hp ' fin=' fn.rtp_op1 ' fout=' fn.rtp_sar1 ...
           ' > /home/chepplew/data/scratch/sar_out.log'];
    command2 = [sartaexe.ha ' fin=' fn.rtp_op2 ' fout=' fn.rtp_sar2 ...
          ' > /home/chepplew/data/scratch/sar_out.log'];
end

tic;   system(command1); toc
tic;   system(command2); toc

% Load in the calcs
[hds1,~,pds1,~] = rtpread(fn.rtp_sar1);
[hds2,~,pds2,~] = rtpread(fn.rtp_sar2);
%
[sfreq, iss] = sort(hds1.vchan);
calc.btcs1 = real(rad2bt(sfreq, pds1.rcalc(iss,:)));
calc.btcs2 = real(rad2bt(sfreq, pds2.rcalc(iss,:)));

% --------------------------------------------------------
%  load the kCARTA calcs. Sort file order first.
% --------------------------------------------------------
% !!!!!!!!! choose PBL (list1) or airslay (list2)  !!!!!!!
disp(['Found ' num2str(length(kc.list1)) ' kCARTA files']);
kc_list = kc.list1;
disp(['Found ' num2str(length(kc.list2)) ' kCARTA files']);
kc_list = kc.list2;

%[ns np] = size(ptemp.rcalc);
% reorder these in profile number order
prfnums=[];
for i=1:length(kc_list)
  fnparts{i} = strsplit(kc_list(i).name, {'_' '.'});
  prfnums(i) = str2num(cell2mat(fnparts{i}(end-1)));
end
[B IB] = sort(prfnums);
kc.rad = []; % zeros(2378,np);
ik = 1;
for i=IB
    x   = load([kc_list(i).folder '/' kc_list(i).name]);
   switch opts.csens
      case {'AIRS_L1C','AIRS_PBL'}
        %krc(:,ik) = x.rKc;
        kc.rad = [kc.rad x.rKc];
        if(i == 1) 
          [kc.freq,kiss] = sort(x.fKc); 
        end
      case {'CRIS_HR','CRIS_PBL'}
        if(isfield(x,'rcris_all'))    kc.rad = [kc.rad x.rcris_all]; end
        if(isfield(x,'hi_rcris_all')) kc.rad = [kc.rad x.hi_rcris_all]; end
        if(i == 1)
          if(isfield(x,'fcris')) i
            [kc.freq, kiss] = sort(x.fcris); end
          if(isfield(x,'hi_fcris')) 
            [kc.freq, kiss] = sort(x.hi_fcris); end
        end
      case {'CHIRP','CHIRP_PBL'}
        kc.rad = [kc.rad x.med_rcris_all];
        if(i == 1) [kc.freq,kiss] = sort(x.med_fcris); end
   end
   ik = ik + 1;
end

calc.rkc1 = kc.rad(kiss,:);
calc.bkc1 = real(rad2bt(kc.freq, kc.rad(kiss,:)));
% REPEAT load kCARTA w/ klist2
calc.rkc2 = kc.rad(kiss,:);
calc.bkc2 = real(rad2bt(kc.freq, kc.rad(kiss,:)));

% ===============================================================================
% Stats and Plots
% ------------------
phome = '/home/chepplew/projects/sarta/prod_2025/figs/';

% choose a subset
usatz = unique(prof.satzen);
% 1. 49 profiles at satzen=0, solzen=150, sfc=black
iiwnt{1} = find(pds1.satzen == usatz(1) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{2} = find(pds1.satzen == usatz(2) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{3} = find(pds1.satzen == usatz(3) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{4} = find(pds1.satzen == usatz(4) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{5} = find(pds1.satzen == usatz(5) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{6} = find(pds1.satzen == usatz(6) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt{7} = find(pds1.satzen == usatz(7) & pds1.solzen == 150 & pds1.emis(1,:) == 1.0);
iiwnt1_7 = [iiwnt{1} iiwnt{2} iiwnt{3} iiwnt{4} iiwnt{5} iiwnt{6} iiwnt{7}];


j = 1;    % satzen angles
figure(1);clf;subplot(211);
  plot(sfreq, mean(calc.btcs1(:,iiwnt{j}),2),'-', sfreq,mean(calc.bkc1(:,iiwnt{j}),2),'-');
  grid on; xlim([600 2800]);legend('sarta','kcarta','location','best')
subplot(212);
  plot(sfreq, mean(calc.btcs1(:,iiwnt{j}),2)-mean(calc.bkc1(:,iiwnt{j}),2),'-');
  hold on;
  plot(sfreq, std(calc.btcs1(:,iiwnt) - calc.bkc1(:,iiwnt),0,2),'-');
  grid on; xlim([600 2800]);legend('mean','std');

%
clf; for j=1:6 
  plot(sfreq, mean(calc.btcs1(:,iiwnt{j}),2)-mean(calc.bkc1(:,iiwnt{j}),2),'-');
  hold on;
end

% ==============================================================
% 6.feb.2025 Sergio expt. trial w/ real obs. 
fn.test1 = ['/asl/s1/sergio/rtp/j1_ccast_hires/allfov//2025/01/08/' ... 
    'cloudy_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp']
% need sym.link for shorter filename
cmd=['ln -s ' fn.test1 ' ip_test_file.rtp'];
system(cmd) 

[head,hatt,prof,patt] = rtpread('ip_test_file.rtp');
[sfreq, iss] = sort(head.vchan);
rad_ham = hamm_app(double(prof.robs1));
bto_a = real(rad2bt(sfreq, rad_ham(iss,:)));
btc_a = real(rad2bt(sfreq, prof.rcalc(iss,:)));
clf; plot(sfreq, mean(bto_a,2),'-')

% run through klayers.airslay
klayersexe.a = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
FIN = 'ip_test_file.rtp';
FOUT = [tempPath '_cris_alay_20250108G099.rtp'];
command_t = [klayersexe.a ' fin=' FIN ' fout=' FOUT ' > /home/chepplew/data/scratch/kla_out.log'];
tic; system(command_t); toc;
% run through CrIS SARTA.airslay
FSAR = [tempPath '_sar_cris_alay_20250108G099.rtp'];
command_s = [sartaexe.ca ' fin=' FOUT ' fout=' FSAR ' > /home/chepplew/data/scratch/sar_out.log'];
tic; system(command_s); toc
[hds,~,pds,~] = rtpread(FSAR);
btcs1a = real(rad2bt(sfreq, pds.rcalc(iss,:)));

% compare with previous CrIS.SARTA.airslay build
sartaexe.cc1='/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev';
command_s = [sartaexe.cc1 ' fin=' FOUT ' fout=' FSAR ' > /home/chepplew/data/scratch/sar_out.log'];
tic; system(command_s); toc
[hds,~,pds,~] = rtpread(FSAR);
btcs1a = real(rad2bt(sfreq, pds.rcalc(iss,:)));
% --------------------------------------
% compare with new CrIS.SARTA.pbl build:
% --------------------------------------
klayersexe.p = '/home/chepplew/gitLib/klayersV205/BinV201/klayers_pbl_wetwater_test';
FIN = 'ip_test_file.rtp';
FOUT = [tempPath '_cris_pblay_20250108G099.rtp'];
command_t = [klayersexe.p ' fin=' FIN ' fout=' FOUT ' > /home/chepplew/data/scratch/kla_out.log'];
tic; system(command_t); toc;
%
FSAR = [tempPath '_sar_cris_pblay_20250108G099.rtp'];
command_s = [sartaexe.cp ' fin=' FOUT ' fout=' FSAR ' > /home/chepplew/data/scratch/sar_out.log'];
tic; system(command_s); toc
[hds,~,pds,~] = rtpread(FSAR);
btcs1p = real(rad2bt(sfreq, pds.rcalc(iss,:)));
% plot(sfreq, mean(bto_a,2)- mean(btcs1p,2),'-')




