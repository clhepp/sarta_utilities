function [calc idx subs] = compare_kc_sar_r49_generic(opts);

% [calc indx] = compare_kc_sar_r49_generic(opts)
%
% Compare calculated TOA radiances from SARTA and kCARTA
%
%  INPUTS: opts.csens [string] 'iasi','cris_hr','cris_lr','airs' 'chirp'
%          opts.prod  [string] production run {'2019','2020'}
%          opts.build [string] run_{build,date} eg 'dec2018','feb2020'
%          opts.regset [string] regression set. {'r49',saf704'}
%
%  OUTPUTS: calc: [structure] with fields: 
%                 bsc (SARTA)  BT calcs
%                 bkc (kCARTA) BT calcs
%
%           idx [structure]: index for different subsets.
%           subs: the subsets
%
%
% prod_2019 includes calcs for AIRS_l1c, CrIS_hr and IASI
%  based off the regular 49 regression profiles, SAF704 &c
%
% CLH. ver 1: restricted to IASI or Cris for regr49 data.
% Mar.2021 CLH. Added CHIRP, AIRS.L1C. Multiple source options
%               and checks.
% Jun 2021 CLH: Add CrIS LR (NSR) option.
% May2023  CLH: changed input params to opts structure. optional additional
%               minor gases to PROF before calc.
% =========================================================================

%{
Test Sets:
Originally output by Sergio to: 
/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/
or /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';

are copied to:
/home/chepplew/data/sarta/prod_2019/generic/
% Fitting kCARTA L2S ODs
od.home = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'
od.path = 'REGR49_400ppm_H2016_Dec2018_AIRS2834/';
od.path = 'REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/';
%}

% ------------------------------------
%                Set Up
% ------------------------------------

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /asl/packages/airs_decon/source          % hamm_app
addpath /home/chepplew/projects/sarta/matlabcode
addpath /home/chepplew/gitLib/rtp_fill_gases     % fill_gasX

disp('Initializing')

% Check input parameters are valid
if( ~all(isfield(opts,{'csens','prod','build','regset'})) 
  error('Please enter 4 option values'); return; 
end

% Check which sensor to compare w/kcarta {'IASI','AIRS','CRIS'}
csens= upper(opts.csens);
allsens = {'AIRS_L1C','CRIS_LR','CRIS_HR','IASI','CHIRP'};
if(~ismember(csens,allsens)) error('Invalid Sensor'); return; end

% check production
prod=lower(opts.prod);
if(~ismember(prod,{'2019','2020','2021','2022','2025'})) 
  error('invalid production date');
  return; 
end
prod_run = ['prod_' prod];

% check build date
build=lower(opts.build);
if(~ismember(build,{'dec2018','feb2020','may2021','jul2022','jan2025a'})) 
  error('invalid build'); 
  return; 
end

% Choose which regression set to compare.
all_regset = {'r49','saf704'};
regset     = opts.regset; % 'r49'; %'saf704';

% check opts.glist
if(~isfield(opts,'glist') )
  disp(['Setting no fill gases']);
  glist = [];
else
  if(~all(ismember(opts.glist, [2 4 5 6 103]) ))
    error('not all valid gases: options are [2 4 5 6]');
    return;
  end
  glist = opts.glist;
  fprintf(1, 'Filling gases: ');fprintf(1,'%4d',glist); fprintf(1,'\n')
end
% glist:

workdir = ['/home/chepplew/projects/sarta/' prod_run '/' lower(csens) '/'];
cd(workdir)

% home for plots
phome = '/home/chepplew/figs_sync/';
phome = ['/home/chepplew/projects/sarta/' prod_run '/' lower(csens) ...
         '/' build '/figs/'];
if(~exist(phome)) mkdir(phome); end

% Set output directory and sarta calc rtp file name
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) ...
         '/' build '/tests/'];
if(~exist(outdr)) mkdir(outdr); end

% Default scan angles, for reference (SARTA limited to < 63.deg).
scang = [0    8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
        65.0428   69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];

% --------------------------------------------------
%   Select the RTP regression profiles to compute 
%      and kCARTA TOA Rads
%   (original set from Sergio - relocated)
% -------------------------------------------------

srcdr  = ['/home/chepplew/data/sarta/' prod_run '/generic/'];

switch prod
  case '2019'
    srcrtp = [srcdr 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp'];
    kpath  = [srcdr 'kcarta/REGR49/400p_1013mb_7angs_unit_seaemis/'];
    %r49_1013_98lev_400p_unitemis_seaemis_7angs_night/14feb2020/'];
    % 2352 = 49prf*8ang*6sfc
    srcrtp = [srcdr 'r49_1013_400p_8angs_sfc_pert.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/surface_pert/'];
%
    srcrtp = [srcdr 'r49_1013_400p_8angs_sfc_pert_v2.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/surface_pert_v2/'];
    kpath  = [srcdr '/kcarta/REGR49/surface_pert_sep20/'];
    srcrtp = [srcdr 'r49_1013_400p_8angs_sfc_pert.rtp'];
%
    srcrtp = [srcdr 'r49_1013_400p_8angs_sfc_pert_815zobs.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/r49_1013_400p_8angs_sfc_pert_815zobs/'];
%
    kpath  = [srcdr 'kcarta/REGR49/400p_1013mb_7angs_unit_seaemis/'];
    srcrtp = [srcdr 'r49_1013_400p_unitemis_seaemis_7angs_night.rtp'];

  case '2022'
    srcrtp = [srcdr 'r49_1013_400p_unitemis_7angs_night.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/r49_1013_400p_unitemis_7angs_night/'];  

    srcrtp = [srcdr 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/r49_1013_98lev_400p_unitemis_seaemis_7angs_night/'];  

    srcrtp = [srcdr 'r49_1100_400p_unitemis_8angs_pert.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/r49_1100_400p_unitemis_8angs_pert/'];

    srcrtp = [srcdr 'r49_1100_400p_seaemis_8angs_pert_v2.rtp'];
    kpath  = [srcdr '/kcarta/REGR49/r49_1100_400p_seaemis_8angs_pert_v2/'];

    srcrtp = [srcdr 'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
    kpath  = [srcdr '/kcarta/SAF704/SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis/'];
   
   case '2025a'
     srcrtp = [srcdr ''];
     kpath =  [srcdr ''];
end

if(~exist(srcrtp)) error('unable to find srcrtp'); return; end;

% ------------------------------------------
% Prepare for processing
% ------------------------------------------
parts1   = strsplit(srcrtp,{'/' '.'});
srcname  = parts1{end-1};

disp('Loading profile data and configuring run')

% Load the RTP profile data into memory.
[head hatt prof patt] = rtpread(srcrtp);
if(isfield(prof,'rcalc')) disp('Removing rcalc');
  prof = rmfield(prof,'rcalc');
end

% rtp_fill_gas() requires something in patt.
patt = struct;
patt = {{'profiles' 'rtime'      'seconds since 0z, 1 Jan 1958'}};

% add HDO depletion to prof.udef(20,:)
if(ismember(103, opts.glist))
  [head, prof, patt] = fill_hdo(head,prof,patt);
end
%
head.pfields = 5;  % robs, model

% Get the subsetting (to be used for plotting)
subs = struct;
subs.ang   = unique(prof.satzen);
subs.uemis = sort(unique(prof.emis(9,:)));  % 0.65 0.75 0.85 0.9451 0.9863 1.00
subs.usolz = sort(unique(prof.solzen));

idx = struct;
jlen = length(prof.satzen)
% switch 686 = 49prf*7ang*2sfc or 2352 = 49prf*8ang*6sfc
switch jlen
  case 343
    for i=1:7
      idx.a{i}  = find(prof.satzen == subs.ang(i));
    end
  case 686
  for i=1:7 
    idx.a{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) == 1); 
    idx.b{i} =  find(prof.satzen == subs.ang(i) & prof.emis(9,:) < 1); 
  end
  subs.names = {'scan angles x 7','2 x surface emissivity: unit and sea'};
  case 2352
  for i=1:8 
    idx.a{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(6)); 
    idx.b{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(5)); 
    idx.c{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(4)); 
    idx.d{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(3)); 
    idx.e{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(2)); 
    idx.f{i} = find(prof.satzen == subs.ang(i) & prof.emis(9,:) == subs.uemis(1)); 
  end
  subs.names = {'8 x scan angles','6 x surface emissivity: a:unity, b:sea .. f:0.5'};
end

% Get the SARTA executable and update the RTP header
% create temporary file space once.
sTemp = mktemp();
fn_rtp1 = [sTemp '_' lower(csens) '_op.rtp'];         % layers file
fn_rtp2 = [sTemp '_' lower(csens) '_sar.rtp'];      % calc file

switch csens
  case 'IASI'
    x      = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq   = x.f_iasi;
    idchan = x.ichan_iasi;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;        % check the kCARTA prediction surface pressure
      szz = size(prof.zobs);
      prof.zobs = 815000.0*ones(1,szz(2));

    outfiles = rtpwrite_12(fn_rtp1,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic_optr';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic_optr_n2o_hno3_so2_nh3_co2';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_hdo_nte';
    %SARTAEXE = '/asl/packages/sartaV108/Bin/sarta_iasi_may09_wcon_nte';  % used in rtp prod.
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_test';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_p2022jul22_dev';
    
  case 'AIRS_L1B'
    hinfo  = hdfinfo('/asl/data/airs/srf/srftables_m140f_withfake_mar08.hdf');
    freq   = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);
    
    rtpwrite(fn_rtp1,head,hatt,prof,patt);
    
    SARTAEXE = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';

  case 'AIRS_L1C'
    hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
    freq = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    %freq    = load('/home/chepplew/myLib/data/airs_f_l1c_lls.txt');
    %idchan  = int32([1:2645])';
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);
    
    rtpwrite(fn_rtp1,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_2834_mar19_basic_optr_co2';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod_v2';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_may19';  
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev'; % or v3
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_xnte_v02';  

  case 'CRIS_LR'
%    x       = load('/home/chepplew/myLib/data/cris_freq_2grd.mat');
    x       = load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    freq    = x.vchan';
    idchan  = x.idchan';
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;
    %np = length(prof.spres); prof.co2ppm = 385*ones(1,np);  
    rtpwrite(fn_rtp1,head,hatt,prof,patt);

    %SARTAEXE = '/asl/packages/sartaV108/BinV201/sarta_crisg4_nov09_wcon_nte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_lrg4_p2021_dev';
        
  case 'CRIS_HR'
    load('/home/chepplew/myLib/data/fcris_hires_4grd.mat','ichan','vchan');
    %freq    = x.vchan;
    %idchan  = x.ichan;
    head.vchan = vchan;
    head.ichan = ichan;
    head.nchan = length(ichan);
    %head.pmax  = 1013;    
    rtpwrite(fn_rtp1,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_basic_optr_co2';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16_nov19A';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16_aug20';  % uses: Data_CrIS_oct16/Coef
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_prod';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev';

  case 'CHIRP'
% use chirp_1702_freq for nonguard channel coefficients, chirp_1702_freq_4grd otherwise.
    x       = load('/home/chepplew/myLib/data/chirp_1702_freq.mat');
      x.ichan = [1:length(x.vchan)]';
%    x       = load('/home/chepplew/myLib/data/chirp_1691_g2.mat');
%    x       = load('/home/chepplew/myLib/data/chirp_1691_freq.mat');
    x       = load('/home/chepplew/myLib/data/chirp_1702_freq_4grd.mat');
    freq    = x.vchan(:);
    idchan  = x.ichan(:);
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);

    rtpwrite(fn_rtp1,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_bbones';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_base_tra_thrm_nte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_p2022jul22_dev';
        
end
%{
junk   = strsplit(SARTAEXE,'/');
fnout  = ['sar_' junk{end} '_' srcname '.rtp'];
sarout = ['outd/' fnout];
%}

% --------------------------------
% Get the SARTA predictions
% --------------------------------
%if(exist('/home/chepplew/logs/sarta/sar_out.log'))
if(~isempty(dir('/home/chepplew/logs/sarta/sar_out*.log')))
  delete('/home/chepplew/logs/sarta/sar_out.log')
end

disp('SARTA running')

switch csens
  case 'IASI'
    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [fn_rtp2 '_1'];  ofn_4 = [fn_rtp2 '_2'];
    tic
      eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sarta/sar_out.log']);
      eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);
    toc
    [hds,~,pds,~] = rtpread_12(ofn_3);
    calc.rsc  = pds.rcalc;
    calc.bsc  = rad2bt(hds.vchan, pds.rcalc);

    % Save the calculations
    %rtpwrite_12(sarout,head,hatt,ptemp,patt);

  case 'AIRS_L1C'
    ifn = fn_rtp1;      ofn = fn_rtp2;  % [tmp '.sar'];
    tic
      %eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
      command=[SARTAEXE ' fin=' ifn ' fout=' ofn ' listp=38,39 > /home/chepplew/logs/sarta/sar_out.log']; 
      system(command);  
    toc
    [hds,~,pds,~] = rtpread(ofn);
    calc.rsc = pds.rcalc;
    calc.bsc = rad2bt(hds.vchan, pds.rcalc);

  case 'CRIS_LR'
    ifn = fn_rtp1;      ofn = fn_rtp2;  % [tmp '.sar'];
    tic
      eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    toc
    [hds,~,pds,~] = rtpread(ofn);
    calc.rsc = pds.rcalc;
    calc.bsc = rad2bt(head.vchan, pds.rcalc);

  case 'CRIS_HR'
    ifn = fn_rtp1;      ofn = fn_rtp2;  % [tmp '.sar'];
    tic
      eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
      %command=[SARTAEXE ' fin=rtpx fout=' sarout ' > /home/chepplew/logs/sarta/sar_out.log']; 
      %system(command);  
    toc
    [hds,~,pds,~] = rtpread(ofn);
    calc.rsc = pds.rcalc;
    calc.bsc = rad2bt(hds.vchan, pds.rcalc);

  case 'CHIRP'
    ifn = fn_rtp1;      ofn = fn_rtp2; % [tmp '.sar'];
    tic
      eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    %command=[SARTAEXE ' fin=rtpx fout=' sarout ' > /home/chepplew/logs/sarta/sar_out.log']; 
    %system(command);  
    toc
    [hds,~,pds,~] = rtpread(ofn);
    calc.rsc = pds.rcalc;
    calc.bsc = rad2bt(hds.vchan, pds.rcalc);

end
disp(['SARTA run took: ']) 


% -------------------------------------------------------------------
% Get the kcarta data
% -------------------------------------------------------------------
klist = dir([kpath 'individual_prof_convolved_kcarta_crisHI_*.mat']);
kc.list = klist;

disp(['Found ' num2str(length(kc.list)) ' kCARTA files']);

%[ns np] = size(ptemp.rcalc);
% reorder these in profile number order
prfnums=[];
for i=1:length(kc.list) 
  fnparts{i} = strsplit(kc.list(i).name, {'_' '.'}); 
  prfnums(i) = str2num(cell2mat(fnparts{i}(end-1)));
end
[B IB] = sort(prfnums);  
kc.rad = []; % zeros(2378,np); 
ik = 1;
for i=IB 
    x   = load([kc.list(i).folder '/' kc.list(i).name]);
    switch csens
      case 'IASI'    
        %krc(:,ik) = x.rKcIasi;
        kc.rad = [kc.rad x.rKcIasi];
        if(i == 1) kc.freq = x.fiasi; end
      case 'CRIS_HR'
        if(isfield(x,'rcris_all'))    kc.rad = [kc.rad x.rcris_all]; end
        if(isfield(x,'hi_rcris_all')) kc.rad = [kc.rad x.hi_rcris_all]; end
	if(i == IB(1)) 
	  if(isfield(x,'fcris')) kc.freq = x.fcris; end
	  if(isfield(x,'hi_fcris')) kc.freq = x.hi_fcris; end
	end
      case 'CRIS_LR'
        if(i == 1) kc.freq = x.lo_fcris; end
	kc.rad = [kc.rad x.lo_rcris_all];

      case 'AIRS_L1C'    
        %krc(:,ik) = x.rKc;
        kc.rad = [kc.rad x.rKc];
        if(i == 1) kc.freq = x.fKc; end
      case 'CHIRP'
        kc.rad = [kc.rad x.med_rcris_all];
	if(i == 1) kc.freq = x.med_fcris; end
    end
    ik = ik + 1;
end

calc.rkc = kc.rad;
calc.bkc = rad2bt(kc.freq, kc.rad);

% Check RTP used to generate kCARTA same as SARTA
kc.rtp = x.use_this_rtp;
disp(kc.rtp)


if(strcmp(csens,'AIRS_L1C'))
  [sfreq, iss]  = sort(hds.vchan);
  [~, isk]  = sort(kc.freq); 
%  [xi  xj] = seq_match(sort(kc.freq), sort(hds.vchan));
  [xi  xj] = seq_match(kc.freq(isk), hds.vchan(iss));
else
  iss = find(hds.vchan);
  isk = iss;
  sfreq = hds.vchan;
end

%{
% -------------- plot checks -------------

freq  = head.vchan(iss);

switch jlen
  case 343        % 7 angles, single sfc emissivity.

    % global mean
    bts_mn = mean(calc.bsc,2);
    btk_mn = mean(calc.bkc,2);
    plot(freq, bts_mn - btk_mn,'-');
      grid on; title('global. 7.angs 1.emis');legend('sar - kc');
      axis([600 2800 -1.5 1.5]);
      % saveas(fh1,[phome 'sar_kc_bias_220719a.png'],'png')

    % satzen means
    clear bts_mn btk_mn; clf
    for i=1:length(idx.a)       % 7-angles
      bts_mn(:,i) = mean(calc.bsc(:,idx.a{i}),2);
      btk_mn(:,i) = mean(calc.bkc(:,idx.a{i}),2);
    end
    hold on; plot(freq, bts_mn - btk_mn,'-');
    jnk = ceil(subs.ang);


  case 686
    % insert NaNs between bands to create gap in plots
    bsc = [NaN(3,686); calc.bsc(4:721,:); NaN(1,686); calc.bsc(722:1377,:); ...
           NaN(1,686); calc.bsc(1378:end,:)];
    bkc = [NaN(3,686); calc.bkc(4:721,:); NaN(1,686); calc.bkc(722:1377,:); ...
            NaN(1,686); calc.bkc(1378:end,:)];
    freq = [sfreq(1:721); NaN; sfreq(722:1594); NaN; sfreq(1595:end)];   % CrIS HR
    freq = [sfreq(1:721); NaN; sfreq(722:1377); NaN; sfreq(1378:end)];   % CHIRP
    [~,iss] = sort(freq); [iss,~] = find(freq); isk = iss;
    % subset e=1 angles 1:6
    iiwnt = [];
    for i=1:6
      iiwnt = [iiwnt idx.a{i}];
    end
    plot(sfreq, mean(calc.bsc(iss,iiwnt),2),'-', sfreq,mean(calc.bkc(isk,iiwnt),2),'-')
  % Each view angle:
    fh1=figure(1);clf;set(gcf,'resize','off');set(fh1,'Position',fh1.Position+[0 0 0 120]);
    h1=subplot(311);hold on;
      for i=1:6 plot(freq, mean(bsc(iss,[idx.a{i}]),2),'-'); end
      grid on;xlim([600 2600]);ylabel('BT (K)');legend('Mean BT','location','north')
      %title([csens ' TOA Mean BT vs view angle'])
    h2=subplot(312);hold on;
    for i=1:6 plot(freq, mean(bsc(iss,[idx.a{i}]),2) - mean(bkc(isk,[idx.a{i}]),2),'-'); end
      grid on;axis([600 2600 -0.5 0.5]);legend('Mean bias','location','north')
      ylabel('BT (K)');
    h3=subplot(313);hold on;
    for i=1:6 plot(freq, std(bsc(iss,[idx.a{i}]) - bkc(isk,[idx.a{i}]),0,2),'-'); end
      grid on;axis([600 2600 -0.1 0.75]);xlabel('wavenumber (cm^{-1})');ylabel('BT (K)');
      legend('std.dev','location','north');
    %% saveas(gcf, [phome 'sar_kc_' lower(csens) '_mean_bias_vs_angle.fig'],'fig')

  case 3920
    % these are the minor gas column pertubation - go to other script.

end

% use first 6 view angles for each surface subset: (a:emis=1, b:emis=sea)
iwnt = [];
for j=1:6
  iwnt = [iwnt idx.b{j}];
end
iwnt = sort(iwnt);

% satzen = 0 and scanang = 0
iwnt = find(prof.satzen(cell2mat(idx.a)) == 0);

plot(head.vchan, mean(calc.bsc(:,iwnt),2) - mean(calc.bkc(:,iwnt),2),'-')

% All angles, e=1, mean and std SARTA - kCARTA
bdiff = calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt);
bbias = mean(calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt),2);
bstdv = std(calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt),0,2);



  radstd = nanstd( calc.rad(:,iwnt) - calc.rkc(:,iwnt),0,2 );
  btbm   = nanmean(calc.bsc(:,iwnt),2);
  mdr    = 1.0*( 1./drdbt(head.vchan,btbm) );
calc.stdv = mdr.*radstd;

fh1=figure(1);set(fh1,'Resize','Off');set(fh1,'Position',fh1.Position+[0 0 0 120]);
  clf(1);h1=subplot(311);
  %plot(head.vchan, mean(calc.bsc(:,idx.a{1}),2),'-', calc.freq,mean(calc.bkc(:,idx.a{1}),2),'-');
  plot(head.vchan(ib), mean(calc.bsc(ib,iwnt),2),'-', ...
       kc.freq(ib),mean(calc.bkc(ib,iwnt),2),'-');
    grid on; legend('SARTA','kCARTA','location','best');xlim([600 2600])
    title([csens ' kC Sar r49 mean'],'interpreter','none');ylabel('BT (K)')
  h2=subplot(312);
    plot(head.vchan(ib), mean(calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt),2),'-');
    grid on;ylabel('dBT K');xlim([600 2600]);legend('SARTA minus kC');
    ylim([-0.55 1.25]);xlabel('wvn cm^{-1}');
  h3=subplot(313);
    plot(head.vchan(ib), calc.stdv,'-','color',[0.5 0.5 0.5]);
    legend('std.dev');xlabel('wavenumber cm^{-1}');ylabel('dBT K');
    grid on;axis([600 2600 -0.05 1.25]);
  %%%saveas(gcf,[phome 'cris_lr_sar_kc_r49_mean_spectrum_sea.fig'],'fig');
    
fh2=figure(2);set(fh2,'Resize','Off');set(fh2,'Position',fh2.Position+[0 0 0 120]);
clf(2);h3=subplot(211);plot(calc.freq, mean(calc.bias(:,cell2mat(idx.a)),2),'-')
   grid on;title([csens ' kCARTA vs SARTA all.angles sea.emis']);legend('mean(kc - sar)');
   h4=subplot(212);plot(calc.freq, std(calc.bias(:,cell2mat(idx.a)),0,2),'-');
   grid on;legend('std(kc - sar)');xlabel('wvn cm^{-1}');ylabel('dBT K');
   
figure(3);clf;plot(calc.freq, mean(calc.bsc(:,idx.a{1}) - calc.bkc(:,idx.a{1}),2),'-')
  grid on; title('SARTA nadir unit.emis minus sea.emis');
  hold on; plot(calc.freq, mean(calc.bkc(:,idx.a{1}) - calc.bkc(:,idx.b{1}),2),'-');
  legend('SARTA','kCARTA')

xLim = [min(calc.freq)-10 max(calc.freq)+10];
fh4=figure(4);clf(4);%set(fh4,'Resize','Off');set(fh4,'Position',fh4.Position+[0 0 0 120]);
  h1=subplot(211);hold on;
    for j=1:7 plot(calc.freq, mean(calc.bias(:,idx.a{j}),2),'-');end
  grid on; legend('mean(sar - kC)','Location','Best');
  h2=subplot(212);hold on;
    for j=1:7 plot(calc.freq, std(calc.bias(:,idx.a{j}),0,2),'-');end
  grid on; legend('std(sar - kc)','9','18','33','45','53','60','Location','Best')
  

  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])

  linkaxes([h3 h4],'x');set(h3,'xticklabel','');
  pp3=get(h3,'position');set(h3,'position',[pp3(1) pp3(2)-pp3(4)*0.1 pp3(3) pp3(4)*1.1])
  pp4=get(h4,'position');set(h4,'position',[pp4(1) pp4(2)+pp4(4)*0.1 pp4(3) pp4(4)*1.1])

  linkaxes([h1 h2 h3],'x');set(h1,'xticklabel','');set(h2,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.05 pp2(3) pp2(4)*1.05])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)-pp2(4)*0.05 pp2(3) pp2(4)*1.05])
  pp3=get(h3,'position');set(h3,'position',[pp3(1) pp3(2)+pp3(4)*0.1 pp3(3) pp3(4)*1.1])

% Clean up band edges (for use with bias and std.dev plots)
grd.chirp   = [1:4 719:723 1371:1378 1702]; 
grd.cris_hr = [1:4 721:722 1594:1595 2235]; 
bnd.cris_hr = [1 721 722 1594 1595 2235];
bnd.chirp   = [1 721 1377 1702];

addpath /asl/matlib/rtptools
  mm_water = mmwater_rtp(hdr,pdr);



%}
%%OPT = -O0 -convert big_endian -extend-source 132 -check all -g -cpp -traceback -fp-stack-check -warn interface
