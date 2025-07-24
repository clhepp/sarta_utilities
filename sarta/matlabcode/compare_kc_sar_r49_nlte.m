% compare_kc_sar_r49_refl_nlte.m
%
% For validating non-LTE SARTA 
% Mar 2021: kCARTA source path updated.
% Aug 2022: update paths, versions and include extended non LTE (xnte)
% Nov 2023: use same opts structure as used for ftc_dev

cd /home/chepplew/projects/sarta/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt, int2bits, mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/chepplew/projects/sarta/matlabcode

% hardwire the kCARTA production run and coefficient calcs (matches existing 
%  directory paths)
csens    = opts.csens;             % 'IASI'; % 'CRIS_LR';
prod_run = ['prod_' opts.prod];
build    = opts.build;             % 'jul2022';
regset   = upper(opts.regset);     %'R49'; %'saf704';

% Choose which regression set to compare.
all_cregr = {'r49','saf704'};

% Enter which sensor to compare w/kcarta {'IASI','AIRS','CRIS'}
allsens = {'AIRS_L1C','CRIS_LR','CRIS_HR','IASI','CHIRP',...
           'AIRS_PBL','CRIS_HR_PBL','IASI_PBL','CHIRP_PBL'};

% home for plots
phome = ['/home/chepplew/projects/sarta/' prod_run '/' lower(csens) ...
         '/' build '/figs/'];
phome = '/home/chepplew/figs_sync/';

% Set output directory and sarta calc rtp file name
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) ...
         '/' build '/tests/'];


% -------------------------------------------------------------------------------
% Location and filename of pertubation test set:
switch build
  case 'dec2018'
    srcdr = ['/home/chepplew/data/sarta/' prod_run '/generic/'];
    fortp = [srcdr 'r49_1013_400p_seaemis_4solz_3satz_.rtp'];
  case 'apr2021'
    srcdr = '/home/chepplew/data/sarta/prod_2021/generic/';
    fortp = [srcdr 'r49_400p_6satzen_19solzen_for_nonlte_v2.rtp'];
  case 'may2021'
    srcdr = '/home/chepplew/data/sarta/prod_2019/generic/';
    fortp = [srcdr 'r49_1013_400p_seaemis_6solz_6satz_.rtp'];
  case 'jul2022'
    srcdr = ['/home/chepplew/data/sarta/' prod_run '/generic/'];
    fortp = [srcdr 'r49_1013_400p_seaemis_19solz_6satz_.rtp'];
  case 'jan2025a'
    srcdr = ['/home/chepplew/data/sarta/' prod_run '/generic/'];
    fortp = [srcdr 'r49_1013_400p_pbl_unitemis_19solz_6satz_.rtp'];

end

% ================================================================================
if(~exist(fortp))

% Original R49 RTP:
  origdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
  srcfn  = [origdr  'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];
   
  [head hatt prof patt] = rtpread(srcfn);
  [nr np] = size(prof.emis);

  % Different view and solar angles combos: 
  % 
  a.satz = [0, 10, 30];
  a.solz = [0, 40, 60, 150];
  %
  b.satz = [0:10:50];
  b.solz = [0:10:80 85 87 90:2:100 105 120];
  %
  c.satz = [0, 10, 20, 30, 40, 50];
  c.solz = [0, 40, 60, 80, 85, 90];
  % New extended nonLTE to 120-deg.
  d.satz = [0, 10, 20, 30, 40, 50];
  d.solz = [0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 87, 90, ...
            92, 94, 96, 98, 100, 105, 120];

  % Now duplicate the 49 profiles for the different angles:
  % SARTA limited to angles < 63-deg
  % *** ONLY do this once to set it up! ***

  % Choose which angles sets
  da = d;

  % Do the satzen and solzen duplication
  h2 = struct;
  p2 = struct;
  nsol = length(da.solz);
  nsat = length(da.satz);

  for ip = 1:48
   [hy,py] = replicate_rtp_headprof(head,prof,ip,nsol*nsat);
   itot = 0;
   for im = 1:nsat
     for in = 1:nsol
      itot = itot+1;
      py.satzen(:,itot) = da.satz(im);
      py.solzen(:,itot) = da.solz(in);
     end
   end
   if ip == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
  end
  % for ease of subsetting and comparing with kCARTA
  idx = struct;
  for i = 1:length(da.solz)
    idx.solz{i} = find(p2.solzen == da.solz(i));
  end
  for i = 1:length(da.satz)
    idx.satz{i} = find(p2.satzen == da.satz(i));
  end
%
% Option to duplicate again with unit emissivity (zero rho)
  unitemis = ones(size(p2.emis));
  zerorho =  zeros(size(p2.rho));
%
  p3 = p2;
  p3.emis  = unitemis;
  p3.rho   = zerorho;
  p2p3 = rtp_cat_prof(p2, p3);

  if(isfield(p2p3,'rcalc')) 
    p2p3 = rmfield(p2p3,'rcalc');
  end
  rtpwrite(fortp, h2, hatt, p2p3, patt);

    clear head prof patt hatt;
end               % ~exist([outgn fortp])

%{
% for ease of subsetting and comparing with kCARTA
idx = struct;
usolz = unique(prof.solzen);
usatz = unique(prof.satzen);
for i = 1:length(usolz)
  idx.solz{i} = find(prof.solzen == usolz(i));
end
for i = 1:length(usatz)
  idx.satz{i} = find(prof.satzen == usatz(i));
end
idx.seaemis  = find(prof.emis(9,:) < 0.995);
idx.oneemis  = find(prof.emis(9,:) == 1);
%}

% ========================================================================= 
% Prepare for the SARTA runs

% create scratch for temporary files:
tempDir = mktemp();

% Use exisiting RTP to do the test:
fnrtp = fortp;
[head hatt prof patt] = rtpread(fnrtp);
if(isfield(prof,'rcalc')) 
  prof = rmfield(prof,'rcalc');
end
prof.spres      = 1013.25 * ones(size(prof.spres));
prof.udef(20,:) = zeros(size(prof.rlat));

% update header for chosen sensor and write rtp for SARTA 
% and select SARTA executable

switch csens
  case {'IASI','IASI_PBL'}
    x      = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq   = x.f_iasi;
    idchan = x.ichan_iasi;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
%    h2.pmax  = 1013;        % check the kCARTA prediction surface pressure
    fn.op_rtp = [tempDir '_' lower(csens) '.rtp'];
    outfiles  = rtpwrite_12(fn.op_rtp, head,hatt,prof,patt);
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_hdo_nte';
    % SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_nov19';
    % SARTAEXE = '/asl/packages/sartaV108/Bin/sarta_iasi_may09_wcon_nte';
    % SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_xnte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jul22_dev_xnte';

  case {'AIRS_L1C','AIRS_PBL'}

    hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf')
    freq  = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    %freq    = load('/home/chepplew/myLib/data/airs_f_l1c_lls.txt');
    %idchan  = int32([1:2645])';
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);
    fn.op_rtp = [tempDir '_' lower(csens) '.rtp'];
    rtpwrite(fn.op_rtp, head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_2834_mar19_basic_optr_tra_nte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_pbl_dev_debug';
    
  case 'CRIS_LR'
    x       = load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    freq    = x.vchan;
    idchan  = x.idchan;
    head.vchan = freq';
    head.ichan = idchan';
    head.nchan = length(idchan);
    %head.pmax  = 1013;    
    fn.op_rtp = [tempDir '_' lower(csens) '.rtp'];
    rtpwrite(fn.op_rtp, head,hatt,prof,patt);
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_basic_optr_co2';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_lrg4_p2021_dev';

  case {'CRIS_HR','CRIS_HR_PBL'}
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    head.pmax  = 1013;    
    fn.op_rtp = [tempDir '_' lower(csens) '.rtp'];
    rtpwrite(fn.op_rtp, head,hatt,prof,patt);
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_basic_optr_co2';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_xnte_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_pbl_dev'

  case {'CHIRP','CHIRP_PBL'}
    x = load('/home/chepplew/myLib/data/chirp_1702_freq.mat');
    head.vchan = single(x.vchan(:));
    head.ichan = int32([1:length(x.vchan)]');
    head.nchan = length(x.vchan);
    fn.op_rtp = [tempDir '_' lower(csens) '.rtp'];
    rtpwrite(fn.op_rtp, head,hatt,prof,patt);
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_base_tra_thrm_nte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_p2022jul22_xnte_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta_f90/bin/chirp_reg_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_pbl_p2025jan25a_dev';
  
end
junk   = strsplit(SARTAEXE,'/');
junk2 = strsplit(fortp,{'/','.'});
fnout  = ['sar_' junk{end} '_' junk2{end-1}];
sarout = ['outd/' fnout];

% ====================================================================
%        Get the SARTA predictions (w/nte)
% ====================================================================
%if(exist('/home/chepplew/logs/sar_out.log'))
if( ~isempty(dir('/home/chepplew/logs/sar_out*.log')) )
  delete('/home/chepplew/logs/sar_out*.log')
end

switch csens
  case {'IASI','IASI_PBL'}
    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [tempDir '_sar_1'];  ofn_4 = [tempDir '.sar_2'];
    tic
    eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ...
          ' > /home/chepplew/logs/sarta/sar_out1.log']);
    eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ...
          ' > /home/chepplew/logs/sarta/sar_out2.log']);
    toc
    cfin = [tempDir '.sar'];
    [~,~,ptemp,~] = rtpread_12(cfin);
    % Save the calculations
    %rtpwrite_12(sarout,head,hatt,ptemp,patt);
    calc.rsc  = ptemp.rcalc;
    calc.bsc  = rad2bt(head.vchan, ptemp.rcalc);

  otherwise     % case {'AIRS','CRIS','CHIRP'}
    ifn = fn.op_rtp;      
    ofn = [tempDir '_sar_' lower(csens) '.rtp'];
    command = [SARTAEXE ' fin=' ifn ' fout=' ofn ...
             ' > /home/chepplew/logs/sarta/sar_out.log']; 
    tic; system(command); toc 
    [~, ~, ptemp, ~] = rtpread(ofn);
    calc.rsc  = ptemp.rcalc;
    calc.bsn  = real(rad2bt(head.vchan, ptemp.rcalc));

end

% Get the SARTA angle subsets (for unity emissivity or sea surface)
iibck   = find(ptemp.emis(1,:) == 1);
iisea   = find(ptemp.emis(1,:) < 0.99);
ssatz=unique(ptemp.satzen);
ssolz=unique(ptemp.solzen);
clear iis;
for i = 1:length(ssatz)
  for j = 1:length(ssolz)
    iis{i}{j} = find(ptemp.satzen(iibck) == ssatz(i) & ...
                     ptemp.solzen(iibck) == ssolz(j));
  end
end

% ===============================================================
%               run SARTA w/o nonLTE 
% ===============================================================
if(exist('/home/chepplew/logs/sar_out.log'))
  delete('/home/chepplew/logs/sar_out.log')
end

switch csens
  case 'IASI'
    SARTAEXE2 = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_nte_off';
    tic
    eval(['! ' SARTAEXE2 ' fin=' ifn_1 ' fout=' ofn_3 ...
          ' > /home/chepplew/logs/sar_out.log']);
    eval(['! ' SARTAEXE2 ' fin=' ifn_2 ' fout=' ofn_4 ...
          ' > /dev/null']);
    toc
    cfin = [tempDir '_sar'];
    [~,~,ptemp2,~] = rtpread_12(cfin);

  otherwise %   case {'AIRS','CRIS_LR','CHIRP','CRIS_HR'}
    ifn = fn.op_rtp;      
    ofn = [tempDir '_sar_' lower(csens) '.rtp'];
    command = [SARTAEXE ' fin=' ifn ' fout=' ofn ...
             ' > /home/chepplew/logs/sarta/sar_out.log']; 
    tic; system(command); toc 
    [~, ~, ptempe, ~] = rtpread(ofn);
    calc.rse  = ptempe.rcalc;
    calc.bse  = real(rad2bt(head.vchan, ptempe.rcalc));

end



% ====================================================================
%                Get the kCARTA TOA rads
% ====================================================================
kc = struct;
% update Notes: for AIRS_L1C (april 2019)
% 
kc.home = '/home/sergio/KCARTA/NONLTE_PRODUCTION/';
kc.home = '/asl/s1/sergio/NLTE_CALCS/';
switch build
  case 'dec2018'
    kc.dir400 = [kc.home 'VT_48PROFILES_120_400ppmv_v120_H16_NLTEH16_Dec2018/'...
                 'Results/CONV_Results/'];
  case 'apr2019'
    kc.dir400 = [kc.home 'VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Apr2019/'...
                 'Results/CONV_Results_Apr23_2019/'];
  case 'jan2020'
    kc.dir400 = [kc.home 'VT_48PROFILES_120_400ppmv_v120_H16_NLTEH16_Dec2018/'...
                 'Results/CONV_Results/'];
    kc.dir385 = [kc.home 'VT_48PROFILES_120_385ppmv_v121_H16_NLTEH16_Jan2020/'...
                 'Results/CONV_Results/'];
  case 'mar2020'
    kc.dir385 = ['/asl/s1/sergio/VT_48PROFILES_120_SRCv1.21_385ppmv_H16_Mar2020/'...
                'CONV_Results/'];
    kc.dir400 = ['/asl/s1/sergio/VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Mar2020/'...
                 'CONV_Results/'];
  case {'apr2021','jul2022'}    % Extended to 120-deg solzen.
    kc.dir385 = [kc.home ...
                 'VT_48PROFILES_120_SRCv1.21_385ppmv_H16_Apr2021_NewNLTEProfiles/'...
                 'CONV_Results/'];
    kc.dir400 = [kc.home ...
                 'VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Apr2021_NewNLTEProfiles/'...
                 'CONV_Results/'];
  case 'may2021'
    kc.dir385 = [kc.home 'VT_48PROFILES_120_385ppmv_v121_H16_NLTEH16_Jan2020/Results/CONV_Results/'];
    kc.dir400 = [kc.home 'VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Jan2020/Results/CONV_Results/'];

end

kc.list3 = dir([kc.dir385 'vt*.mat']);
kc.list4 = dir([kc.dir400 'vt*.mat']);
kc.list  = kc.list4;

% reorder these in profile number order
clear prfnums;
for i=1:length(kc.list) 
  fnparts{i} = strsplit(kc.list(i).name, {'vt' '.'}); 
  prfnums(i) = str2num(cell2mat(fnparts{i}(2)));
end
[IB,~] = sort(prfnums);  
% IBX = IB([1:9 12:48]);   % NewNLTE file:10,11 is corrupt
kc.re = []; kc.rn = []; kc.solz = []; kc.satz = [];
for i=IB([1:9 12:48]) 
  x   = load([kc.list(i).folder '/' kc.list(i).name]);
  switch csens
    case 'IASI'    
      kc.rn = [kc.rn x.nlte_iasiB];
      kc.re = [kc.re x.lte_iasiB];
      kc.freq  = x.fiasi;
    case 'CRIS_LR' 
      kc.rn = [kc.rn x.lo_nlte_crisB];
      kc.re = [kc.re x.lo_lte_crisB]; 
      kc.freq  = x.lo_fcris;
    case 'CRIS_HR' 
      kc.rn = [kc.rn x.hi_nlte_crisB];
      kc.re = [kc.re x.hi_lte_crisB]; 
      kc.freq  = x.hi_fcris;
    case 'AIRS_L1C'    
      kc.rn = [kc.rn x.nlte_airsB];
      kc.re = [kc.re x.lte.airsB]; 
      kc.freq  = x.fairs;
    case 'CHIRP'    
      kc.rn = [kc.rn x.med_nlte_crisB]; 
      kc.re = [kc.re x.med_lte_crisB];
      kc.freq  = x.med_fcris;
  end
  kc.solz  = [kc.solz x.solar];
  kc.satz  = [kc.satz x.viewer];
  fprintf(1,'.');
end

calc.bkn = real(rad2bt(kc.freq, kc.rn));
calc.bke = real(rad2bt(kc.freq, kc.re));

% Get kCARTA angle subsets
ksatz=unique(kc.satz);
ksolz=unique(kc.solz);
clear iik;
for i = 1:length(ksatz)
  for j = 1:length(ksolz)
    iik{i}{j} = find(kc.satz == ksatz(i) & kc.solz == ksolz(j));
  end
end

% ---------------------------------------------------------------
%  Alternative kCARTA predict for absolute night calc
% ----------------------------------------------------------------
kpath2='/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Dec2018_AIRS2834/';
prfnums2=[]; clear fnparts;
%if(strcmp(cregr,'r49'))
  fnlst2 = dir([kpath2 'RAD1013_seaemis_Dec13_2019/convolved_kcarta_RAD1013_*_radiances.mat']);
% reorder these in profile number order
  for i=1:length(fnlst2) 
    fnparts{i}  = strsplit(fnlst2(i).name, '_');
    prfnums2(i) = str2num(cell2mat(fnparts{i}(4)));
  end
  [B IB] = sort(prfnums2);
  krc2 = [];
  for i=IB
    x   = load([fnlst2(i).folder,'/' fnlst2(i).name]);
    if(strcmp(csens,'IASI'))    krc2 = [krc2 x.riasi_all]; end % 49 profs x 8 angles
    if(strcmp(csens,'CRIS_HR')) krc2 = [krc2 x.rcris_all]; end
    if(strcmp(csens,'AIRS'))    krc2 = [krc2 x.rairs_all]; end
  end
  kc2_freq = x.fiasi;

krc2_nadir = krc2(:,1:8:392);
bkc2_nadir = real(rad2bt(kc2_freq, krc2_nadir));


% ============================================================= 
%                 PLOTTING SECTION
% ---------------------------------------------------------
%
addpath /asl/packages/airs_decon/source
if(contains(csens,'AIRS'))
  [sfreq,  iia] = sort(head.vchan);
  [kfreq,  iik] = sort(kc.freq);
  [ia, ib] = seq_match(sfreq, kfreq);
else
  iia = find(head.vchan);
  iib = find(kc.freq);
end

junk2 = round(1000*head.vchan)./1000;
junk1 = round(1000*kc.freq)./1000;
[aaa, bbb] = intersect(sort(junk1), sort(junk2));
% KCARTA variation with solzen, satzen=0
figure(2);clf;hold on;
  plot(kc.freq(iib), mean(calc.bkn(iib,iik{1}{1}),2),'-'); 
  plot(kc.freq(iib), mean(calc.bkn(iib,iik{1}{3}),2),'-');
  plot(kc.freq(iib), mean(calc.bke(iib,iik{1}{1}),2),'-');
  grid on; xlabel('wvn cm^{-1}');ylabel('BT (K)');
  xlim([2110 2500]);legend('Noon','Sunset','Location','Best');
  title([csens ' kCARTA nonLTE mean regr49'],'interpreter','none')

% SARTA variation with solzen, satzen=0
figure(3);clf;hold on;
  for j=[1:2:12 16 19]
    plot(sfreq, mean(calc.bsn(iia,iis{1}{j}),2),'-'); 
  end
  plot(sfreq, mean(calc.bsn(iia,iis{1}{3}),2),'-');
  plot(sfreq, mean(calc.bse(iia,iis{1}{1}),2),'-');
  grid on; xlabel('wvn cm^{-1}');ylabel('BT (K)');
  xlim([2110 2500]);legend('Noon','Night','Location','Best');
  title([csens ' SARTA nonLTE mean regr49']),'interpreter','none'

% SARTA vs kCARTA Day-Night Difference
figure(2);clf;hold on;
  plot(head.vchan(iia), mean(calc.bsn(iia,iis{1}{1}),2) - ...
                        mean(calc.bsn(iia,iis{1}{4}),2),'-');
  plot(kc.freq(iib), mean(calc.bkn(iib,iik{1}{1}),2) - ...
                     mean(calc.bkn(iib,iik{1}{6}),2),'-');
  grid on; xlabel('wvn cm^{-1}');ylabel('BT (K)');
  xlim([2000 2780]);legend('Noon - Night');
  title([csens ' SARTA nonLTE'],'interpreter','none')

% SARTA vs kCARTA nonLTE - LTE Difference
figure(4);clf;hold on;
  plot(head.vchan(iia), mean(calc.bsn(iia,iis{1}{1}),2) - ...
                        mean(calc.bse(iia,iis{1}{1}),2),'-');
  plot(kc.freq(iib),    mean(calc.bkn(iib,iik{1}{1}),2) - ...
                        mean(calc.bke(iib,iik{1}{1}),2),'-');
  grid on; xlabel('wvn cm^{-1}');ylabel('BT (K)');
  xlim([2000 2780]);legend('SARTA','kCARTA');
  title([csens ' SARTA nonLTE - LTE'],'interpreter','none')



% ------- Variation with solzen, select channels --------
% rms a group of channels
freq = kc.freq;
ichs = find(freq > 2351,4);

for k = 1:length(iik)
  for j = 1:length(iik{1})
    kc.rms_n(k,j) = rms(calc.bkn(ichs, iik{k}{j}),[1 2]);
    kc.rms_e(k,j) = rms(calc.bke(ichs, iik{k}{j}),[1 2]);
  end
end
%
for k = 1:length(idx.satz)
  for j = 1:length(idx.solz)
    iiwnt = intersect(idx.seaemis, intersect(idx.satz{k}, idx.solz{j}));
    sa.rms_n2(k,j) = rms(calc.bsn2(ichs, iiwnt),[1 2]);
    sa.rms_e(k,j) = rms(calc.bse(ichs, iiwnt),[1 2]);
  end
end

figure;clf;subplot(211);hold on;grid on;ylim([235 262])
  plot(ksolz,kc.rms_n,'m.-', ksolz,kc.rms_e,'c.-')
    subplot(212);hold on;grid on;ylim([235 262])
  plot(usolz,sa.rms_n,'m.-', usolz,sa.rms_e,'c.-')
  legend('kC nlte','kC eq','sar nlte','sar eq')
  xlabel('solzen (deg)');ylabel('BT (K)')
  title('CrIS.LR BT vs Solar Zen 2352 cm-1')

for j=1:7
  calc.bias(j,:) = mean(calc.bkc(:,j:7:343) - calc.bsc(:,j:7:343),2);
  calc.stdv(j,:) = std(calc.bkc(:,j:7:343) - calc.bsc(:,j:7:343),0,2);
end

% ============================================================
%  subset for each solzen for bin averaging and std.dev 
usolz = unique(ptemp.solzen);
usolz_mn = (usolz(2:end) + usolz(1:end-1))/2.0;
clear iisz btbias btstd
for i=1:length(usolz)-1
  iisz{i} = find(ptemp.solzen > usolz(i) & ptemp.solzen <= usolz(i+1));
end
for i=1:length(usolz)-1
  btbias(i) = mean(calc.bsn(iss(ich),iisz{i}) - calc.bse(iss(ich),iisz{i}));
  btstd(i)  = std(calc.bsn(iss(ich),iisz{i}) - calc.bse(iss(ich),iisz{i}),0,2);
end
 
 
  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])

%From Sergio:
%"remember kCARTA NLTE calcs literally uses kCARTA LBL  computed using ancient linemixing
%In other words you are fitting (KCARTA LBL NLTE-KCARTA LBL LTE) for delta(NLTE rad)
%while you are fitting KCARTA LBLRTM LTE for SARTA LTE, and then adding in the fitted delta(NLTErad)
%better just to compare SARTA (NLTE-LTE) against KCARTA LBL (NLTE-LTE)"


