% compare_kc_sar_r49_generic_perturb
%
% validate minor gas perturbation on TOA radiance vs kCARTA
% INPUT:
%  opts structure with fields:
%       csens  [string]
%       prod   [string]
%       build  [string]
%       regset [string]
%       vers   [integer] : the perturbation version (0,1,2)
%
% For perturbations ordering see: ~/projects/sarta/matlabcode/make_r49_perturbed_profs.m 
%      indx.
% CLH
% Jul.2024 CLH: more generic, with opts structure for inputs 
% Jan 2025 CLH: add new paths for PBL variants

cd /home/chepplew/projects/sarta/

addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/motteler/shome/airs_decon/source   % seq_match

%
if(~all(isfield(opts,{'csens','build','prod','regset','vers'})))
  error('Need all FIVE options')
  return;
end

%
csens    = upper(opts.csens); % 'CRIS_HR'; % 'IASI'; % 'CRIS_LR'; % 'CHIRP';
prod_run = ['prod_' opts.prod];
build    = lower(opts.build);   % 'jul2022'; % 'may2021';
regset   = lower(opts.regset);

% Initialize
allsens = {'AIRS_L1C','CRIS_LR','CRIS_HR','IASI','CHIRP',...
           'AIRS_PBL','CRIS_HR_PBL','IASI_PBL','CHIRP_PBL'};
if(~ismember(csens,allsens)) error('Invalid Sensor');
  return; end

% Set output directory and sarta calc rtp file name
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) ...
         '/' build '/tests/'];

% Set output directgory for plots:
phome = ['/home/chepplew/projects/sarta/' prod_run '/' lower(opts.csens) ...
         '/' opts.build '/figs/'];

% since rtpwrite can't handle long paths - ensure symlink outd is correct
if(exist('outd')) disp(['deleting existing symlink for data outd'])
  delete('outd');
end
if(~exist('outd'))
  disp(['creating symlink for data']);
  command = ['ln -s ' outdr ' outd'];
  system(command);
end

% reset results structure
calc = struct;

% original location from Sergio:
dpath='/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/';
% and /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/ ...
%  r49_1100_400p_unitemis_nadir_pert.rtp
%  r49_1100_400p_unitemis_8angs_pert.rtp
% Copied location to ~chepplew/data/sarta/prod_2019/generic/:

% -----------------------------------------------------------
% !!!!!!!!!!! switch for version of perturbations !!!!!!!!!!!
swtich opts.vers
  case 0
% -------  nadir set: -------
% Original profile RTP file used to generate TOA radiances.
    fnrtp = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
             'r49_1100_400p_unitemis_nadir_pert.rtp'];
    kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
             '/kcarta/REGR49/r49_1100_400p_unitemis_8angs_pert/'];

% This matrix gives the order of profiles after duplication and perturbation:
    pind =  {{'WV',   1,  [1:9:433]},...
            {'CO2',  2,  [2:9:434]},...
            {'O3',   3,  [3:9:435]},...
            {'N2O',  4,  [4:9:436]},...
            {'CO',   5,  [5:9:437]},...
            {'CH4',  6,  [6:9:438]},...
            {'SO2',  9,  [7:9:439]},...
            {'HN3', 11,  [8:9:440]},...
            {'NHO3',12,  [9:9:441]},...
            {'unp',  0,  [442:490]}};

% -------------------------------------------------
% Algorithm for indexing each 8 angle perturbation (49*8*10=3920)
% -------------------------------------------------
% Original profile RTP file used to generate TOA radiances.
fnrtp = ['/home/chepplew/data/sarta/prod_2019/generic/' ...
         'r49_1100_400p_unitemis_8angs_pert.rtp'];
fnrtp = ['/home/chepplew/data/sarta/prod_2020/generic/' ...
         'r49_1100_400p_unitemis_8angs_gas_pert_v1.rtp'];
fnrtp = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
         'r49_1100_400p_unitemis_8angs_pert.rtp'];
kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
         '/kcarta/REGR49/r49_1100_400p_unitemis_8angs_pert/'];

  case 1
% -------------------------------
% Version 1: simple 10% all gases
% -------------------------------
    if(contains(csens,'PBL'))
      fnrtp = ['/home/chepplew/data/sarta/prod_2025/generic/' ...
              'r49_1013_400p_pbl_unitemis_8angs_gas_pert_v1.rtp'];
      kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
            'kcarta/r49_1013_400p_pbl_unitemis_8angs_gas_pert_v1/'];
    else
      fnrtp = [''];
      kpath = [''];
    end
    indx = struct;
    for ii=1:49 
      for kk=1:8 
        jj=8*(ii-1)+kk;
        indx.wv(jj)   = kk    + (ii-1)*72;
        indx.co2(jj)  = 8+kk  + (ii-1)*72;
        indx.o3(jj)   = 16+kk + (ii-1)*72;
        indx.n2o(jj)  = 24+kk + (ii-1)*72;
        indx.co(jj)   = 32+kk + (ii-1)*72;
        indx.ch4(jj)  = 40+kk + (ii-1)*72;
        indx.so2(jj)  = 48+kk + (ii-1)*72;
        indx.nh3(jj)  = 56+kk + (ii-1)*72;
        indx.hno3(jj) = 64+kk + (ii-1)*72;
        %disp([num2str(jj) ': ' num2str(indx.so2(jj))]);
      end 
    end
    indx.unp = [3529:3920];

  case 2
% -----------------------
% version 2-Mar-2019 LLS 
% -----------------------
% Original profile RTP file used to generate TOA radiances.
    fnrtp = ['/home/chepplew/data/sarta/prod_2019/generic/' ...
             'r49_1100_400p_seaemis_8angs_pert_v2.rtp'];
    kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
             '/kcarta/REGR49/r49_1100_400p_seaemis_8angs_pert_v2/'];
%
    fnrtp = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
             'r49_1100_400p_unitemis_8angs_pert.rtp'];
    kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
             '/kcarta/REGR49/r49_1100_400p_unitemis_8angs_pert/'];

if(contains(csens,'PBL') & opts.build == 'jan2025a'))
   fnrtp = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
            'r49_1013_400p_pbl_unitemis_8angs_gas_pert_v1.rtp'];
   kpath = ['/home/chepplew/data/sarta/' prod_run '/generic/' ...
            'kcarta/r49_1013_400p_pbl_unitemis_8angs_gas_pert_v2/'];
end

    indx = struct;
    for ii=1:49 
      for kk=1:8 
        jj=8*(ii-1)+kk;
        indx.wv(jj)   = kk    + (ii-1)*72;
        indx.co2(jj)  = 8+kk  + (ii-1)*72;
        indx.o3(jj)   = 16+kk + (ii-1)*72;
        indx.n2o(jj)  = 24+kk + (ii-1)*72;
        indx.co(jj)   = 32+kk + (ii-1)*72;     %  CO + 10%
        indx.cox(jj)  = 40+kk + (ii-1)*72;     %  CO + 30%
        indx.so2(jj)  = 48+kk + (ii-1)*72;
        indx.nh3(jj)  = 56+kk + (ii-1)*72;
        indx.hno3(jj) = 64+kk + (ii-1)*72;
        %disp([num2str(jj) ': ' num2str(indx.so2(jj))]);
        
      end 
    end
    indx.unp = [3529:3920];
end
% !!!!!!!!!!! end switch for version of perturbations !!!!!!!!!!!
% -----------------------------------------------------------

% -----------------------------------------------------------
% Get and modify header for requested sensor & select SARTA executable:
[head hatt prof patt] = rtpread(fnrtp);
if(isfield(prof,'rcalc'))
  prof = rmfield(prof,'rcalc');
end
% set sfc pressure 1013 mb
%[~,ns] = size(prof.spres)
%prof.spres = 1013.250* ones(1,ns);
head.pmax = 1013.0;
prof.udef(20,:) = zeros(size(prof.rlat));

% create scratch/temporary folder (once)
tempDir = mktemp();

switch csens
  case 'CRIS_LR'
    x       = load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    freq    = x.vchan;
    idchan  = x.idchan;
    head.vchan = freq';
    head.ichan = idchan';
    head.nchan = length(idchan);
    rtpwrite(tmp,head,hatt,prof,patt);

    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_lrg4_p2021_dev';

  case {'CRIS_HR','CRIS_HR_PBL'}
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    fn.op_csens = [tempDir '_cris_hr_pbl_pert_2235_op.rtp'];
    rtpwrite(fn.op_csens, head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_basic_optr_co2_n2o';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_pbl_dev';

  case {'AIRS_L1C','AIRS_PBL'}
    
    hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf')
    freq = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    %freq    = load('/home/chepplew/myLib/data/airs_f_l1c_lls.txt');
    %idchan  = int32([1:2645])';
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);
    fn.op_csens = [tempDir '_airs_pbl_pert_2834_op.rtp'];
    rtpwrite(fn.op_csens, head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_2834_mar19_basic_optr_co2';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_pbl_dev_debug';
  
  case {'IASI','IASI_PBL'}
    x = load('/home/chepplew/gitLib/ftc_dev/chanLists/iasi_8461_list');
    head.vchan = single(x(:,2));
    head.ichan = int32(x(:,1));
    head.nchan = length(x);    
    %fnrtpx = ['outd/r49_1013m_400p_e1_nadir_pert_8461.rtp'];
    fnrtpx = ['outd/r49_1013m_400p_e1_8ang_pert_8461.rtp'];
    fn.op_csens = [tempDir '_iasi_pbl_pert_8461_op.rtp'];
    outfiles = rtpwrite_12(fn.op_csens, head,hatt,prof,patt);
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic_optr_n2o_hno3_so2_nh3_co2';
    %SARTAEXE = '/asl/packages/sartaV108/BinV201/sarta_iasi_may09_wcon_nte';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jul22_dev';

  case {'CHIRP','CHIRP_PBL'}
    x  = load('/home/chepplew/myLib/data/chirp_1702_freq_4grd.mat');  % 1702 w/grd
    head.vchan = single(x.vchan(:));   % single(x.wnum(:));
    head.ichan = int32(x.ichan(:));    % int32([1:length(x.wnum)]');
%    x  = load('/home/motteler/shome/chirp_test/chirp_wnum.mat');    % 1679 grid
%    head.vchan = single(x.wnum(:));
%    head.ichan = int32([1:length(x.wnum)]');
    head.nchan = length(head.ichan);    
    fn.op_csens = [tempDir '_chirp_pbl_pert_1679_op.rtp'];
    rtpwrite(fn.op_csens, head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_base_so2_nh3_hno3_n2o';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_p2022jul22_xnte_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_pbl_p2025jan25a_dev';
end

% ------------------------------------------------
% Get the SARTA TOA RAD predicts:
% ------------------------------------------------

junk   = strsplit(SARTAEXE,{'/','_'});
fnout  = ['_sar_' junk{end} '_8apg_v' num2str(opts.vers) '.rtp'];
fn.sar_csens = [tempDir '_sar_' lower(csens) '.rtp'];

% Need short file names
if(exist('rtpx')) disp(['removing existing sym link']);
  delete('rtpx');
end
command=['ln -s ' fnrtpx ' rtpx'];
system(command);

% clear out old SARTA log files:
if(exist('/home/chepplew/logs/sar_out.log')) 
  disp('Deleting previous SARTA log file');
  delete('/home/chepplew/logs/sar_out.log')
end

switch csens
  case 'IASI'
    ifn_1 = outfiles{1};      ifn_2 = outfiles{2};
    ofn_3 = [tmp '.sar_1'];   ofn_4 = [tmp '.sar_2'];

    %eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' listp=1,2,3,4,5,6,7,8,9']);  % ' > /dev/null']);
    eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sar_out3.log']);
    eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);

    cfin = [tmp '.sar'];
    [~,~,ptemp,~] = rtpread_12(cfin);
    calc.bsc = rad2bt(head.vchan, ptemp.rcalc);

  otherwise    % all other sensors

    ifn = fn.op_csens;      ofn = fn.sar_csens;
    %command=[SARTAEXE ' fin=' fnrtpx ' fout=sar_out.rtp >/dev/null 2>&1']
    %command=[SARTAEXE ' fin=' fnrtpx ' fout=sar_out.rtp > /home/chepplew/logs/sarta/sar_out.log']
    command=[SARTAEXE ' fin=' ifn ' fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log'];
    tic;  system(command); toc
    [hds, ~, pds, ~] = rtpread(ofn);
    [sfreq, iss] = sort(hds.vchan);
    calc.rsc = pds.rcalc(iss,:);
    calc.bsc = real(rad2bt(sfreq, pds.rcalc(iss,:)));

end
% Save the calculations
%rtpwrite_12(sarout,head,hatt,ptemp,patt);

% -----------------------------------------------------------
% Get the kCARTA TOA RAD predicts: NB get correct version!!!
% -----------------------------------------------------------
klist = dir([kpath 'individual_prof_convolved_kcarta_crisHI_*']);
%
klist = dir([kpath 'convolved_kcarta_*.mat']);
%
klist = dir([kpath 'individual_prof_convolved_kcarta_crisHI_crisMED*.mat']);

klist = dir([kpath 'individual_prof_convolved_kcarta_airs_iasi_crisHI_crisMED*.mat']);

% reorder these in profile number order
clear prfnums fnparts;
for i=1:length(klist)
  fnparts{i} = strsplit(klist(i).name, {'_' '.'});
  % prfnums(i) = str2num(cell2mat(fnparts{i}(7)));
  prfnums(i) = str2num(cell2mat(fnparts{i}(9)));
end
[IB IC] = sort(prfnums);

[ns, np] = size(pds.rcalc);
calc.rkc = []; % zeros(ns, np);
calc.bkc = []; % zeros(ns, np);

% NB: CrIS field names changed after adding CHIRP
k = 1;
for i = IC %% [1:length(klist)]
 %disp([klist(i).name])
  switch csens
  case 'CRIS_LR'
    x = load([klist(i).folder '/' klist(i).name]);
    calc.rkc(:,k) = x.lo_rcris_all;
    if(k==1) kc_freq = x.lo_fcris; end
  case {'CRIS_HR','CRIS_HR_PBL'}
    x = load([klist(i).folder '/' klist(i).name]);
    if(isfield(x,'fcris'))    calc.rkc(:,i) = x.rcris_all; end
    if(isfield(x,'hi_fcris')) calc.rkc(:,k) = x.hi_rcris_all; end
    if(k==1 & isfield(x,'hi_fcris')) kc_freq = x.hi_fcris; end
  case {'AIRS_L1C','AIRS_PBL'}
    x = load([klist(i).folder '/' klist(i).name]);
    calc.rkc(:,k) = x.rKc;
    if(k==1) kc_freq = x.fKc; end
  case {'IASI','IASI_PBL'}
    x = load([klist(i).folder '/' klist(i).name]);
    calc.rkc(:,k) = x.rKcIasi;
    if(k==1) kc_freq = x.fiasi; end
  case {'CHIRP','CHIRP_PBL'}
    x = load([klist(i).folder '/' klist(i).name]);
    calc.rkc(:,k) = x.med_rcris_all;
    if(k==1) kc_freq = x.med_fcris; end
  end    
  k = k + 1;
  if(~mod(i,100)) fprintf(1,'.'); end
end

% care for mis-matched spectral grids
[kfreq, ikk] = sort(kc_freq);
calc.bkc  = real(rad2bt(kfreq, calc.rkc(ikk,:)));

[ix iy] = seq_match(sfreq, kfreq);
calc.bbias = calc.bsc(ix,:) - calc.bkc(iy,:);
bfreq = sfreq(ix);

%{
% -----------------------------------------------
%               PLOTTING SECTION 
% -----------------------------------------------
if(strcmp(csens,'AIRS_L1C'))
  [bx,ix] = sort(calc.freq);
else
  bx=':'; ix=find(calc.freq);
end

% ------------ 49 set pertubs at nadir: --------------
izo = pind{10}{3};     % unperturbed indexes
ixo = pind{2}{3};      % CO2
ixo = pind{3}{3};      % O3
ixo = pind{4}{3};      % N2O
ixo = pind{5}{3};      % CO
ixo = pind{6}{3};      % CH4
ixo = pind{7}{3};      % SO2
ixo = pind{8}{3};      % 8 NH3
ixo = pind{9}{3};      % HNO3
bias.kc  = calc.bkc(:,ixo) - calc.bkc(:,izo);
bias.sar = calc.bsc(:,ixo) - calc.bsc(:,izo);
   plot(calc.freq, mean(bias.kc,2),'-', calc.freq, mean(bias.sar,2),'-')
figure(2);clf;plot(calc.freq, mean(calc.bkc(:,izo),2),'-');
      hold on;plot(calc.freq, mean(calc.bsc(:,izo),2),'-');
figure(2);clf;plot(calc.freq, mean(calc.bkc(:,izo),2)-mean(calc.bsc(:,izo),2),'-');


% -------------- 3920 set (49*8*10) ----------------------
% H2O
fign = 1;
for j=1:8 
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.wv(j:8:392); % nadir
  bias(j).kc_h2o  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_h2o = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_h2o  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_h2o = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
figure(fign);clf;hold on; grid on;
   plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
   plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
figure(fign);clf;hold on; grid on; xlim([600 2780]);
   plot(kfreq, mean(bias(j).kc_h2o,2),'-');
   plot(sfreq, mean(bias(j).sar_h2o,2),'-'); 
   plot(kfreq, stdv(j).kc_h2o,'-',sfreq,stdv(j).sar_h2o,'-');
   plot(sfreq(ix), mean(bias(j).sar_h2o(ix,:),2) - mean(bias(j).kc_h2o(iy,:),2),'-');
   %plot(calc.freq(ix), mean(bias.sar_h2o(ix,:) - bias.kc_h2o(ix,:),2),'-');
   %plot(calc.freq(ix), std(bias.kc_h2o(ix) - bias.sar_h2o(ix),0,2),'-')
   th1=title([csens ' Mean,std,bias d(BT) 10% H2O.pert']);set(th1,'Interpreter','none'); 
   legend('kCARTA mean','SARTA mean','kCARTA std','SARTA std','mean bias');
   xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
   %saveas(gcf, [phome 'kc_sar_pert_10pc_h2o_0deg.fig'],'fig');

% pCO2
fign = fign+1;
for j=1:7    % j=1:nadir. j=2: 8.8-deg
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.co2(j:8:392); 
  bias(j).kc_co2  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_co2 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_co2  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_co2 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;hold on;
    plot(kfreq, mean(calc.bkc(:,xzo{j}),2),'-')
    plot(sfreq, mean(calc.bsc(:,xzo{j}),2),'-')
 figure(fign);clf;hold on;grid on;xlim([640 820])
    plot(kfreq(ix), mean(bias(j).kc_co2(ix),2),'-');
    plot(sfreq(ix), mean(bias(j).sar_co2(ix),2),'-'); 
    plot(kfreq(ix), stdv(j).kc_co2(ix),'-', sfreq(ix), stdv(j).sar_co2(ix),'-');
    plot(sfreq(ix), mean(bias(j).sar_co2(ix,:),2) - mean(bias(j).kc_co2(ix,:),2),'-');
    %plot(calc.freq, std(bias.kc_co2 - bias.sar_co2,0,2),'-')
    th1=title([csens ' Mean,std,bias d(BT) from 10% CO2.pert']);set(th1,'Interpreter','none'); 
    legend('kCARTA mean','SARTA mean','kCARTA std','SARTA std','std kc-sar');
    xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
    xlim([600 2600])
    %saveas(gcf, [phome 'kc_sar_pert_10pc_co2_0deg.fig'],'fig');

% Ozone
fign = fign+1;
for j=1:7    % j=1:nadir. j=2: 8.8-deg
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.o3(j:8:392); 
  bias(j).kc_o3  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_o3 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_o3  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_o3 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;hold on;
    plot(kfreq, mean(calc.bkc(:,xzo{j}),2),'-')
    plot(sfreq, mean(calc.bsc(:,xzo{j}),2),'-')
 figure(fign);clf;hold on;grid on;  xlim([640 1200])
 plot(kfreq(iy), mean(bias(j).kc_o3(iy),2),'-');
    plot(sfreq(ix), mean(bias(j).sar_o3(ix),2),'-'); 
    plot(kfreq(iy), stdv(j).kc_o3(iy),'-', sfreq(ix), stdv(j).sar_o3(ix),'-');
    %plot(sfreq, std(bias.kc_co2 - bias.sar_o3,0,2),'-')
    plot(sfreq(ix), mean(bias(j).sar_o3(ix,:),2) - mean(bias(j).kc_o3(iy,:),2),'-');
    th1=title([csens ' Mean,std,bias d(BT) 10% O3.pert']);set(th1,'Interpreter','none'); 
    legend('kCARTA mean','SARTA mean','kCARTA std','SARTA std','std kc-sar');
    xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)');
    %saveas(gcf, [phome 'kc_sar_pert_10pc_o3_0deg.fig'],'fig');

% N2O
fign = fign + 1;
for j=1:7    % j=1:nadir. j=2: 8.8-deg
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.n2o(j:8:392); 
  bias(j).kc_n2o  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_n2o = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_n2o  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_n2o = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j = 1;
 figure(fign);clf;plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
  figure(fign);clf; hold on; grid on; xlim([1150 2600]);
    plot(kfreq(iy), mean(bias(j).kc_n2o(iy,:),2),'-');
    plot(sfreq(ix), mean(bias(j).sar_n2o(ix,:),2),'-'); 
    plot(kfreq, stdv(j).kc_n2o,'-',sfreq,stdv(j).sar_n2o,'-');
    %plot(kfreq, stdv.kc_n2o,'-', sfreq, stdv.sar_n2o,'-');
    plot(sfreq(ix), mean(bias(j).sar_n2o(ix,:),2) - mean(bias(j).kc_n2o(iy,:),2),'-');
    plot(kfreq(iy), std(bias(j).kc_n2o(iy) - bias(j).sar_n2o(ix),0,2),'-')
    title([csens ' Mean,std,bias d(BT) from 10% n2o.pert'],'Interpreter','none'); 
    legend('kCARTA mean','SARTA mean','kCARTA std.dev','SAR std.dev');
    xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
    %saveas(gcf, [phome 'kc_sar_pert_10pc_n2o_0deg.fig'],'fig');

% CH4
fign = fign + 1;
for j = 1:7
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.ch4(j:8:392);
  bias(j).kc_ch4  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_ch4 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_ch4  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_ch4 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j = 1;
 figure(fign);clf;plot(calc.freq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(calc.freq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
      legend('kCARTA','SARTA');
fh=figure(fign);clf;hold on; grid on;  xlim([1180 1650])
%set(gcf,'Resize','Off');set(fh,'Position',fh.Position+[0 0 0 120]);
   plot(calc.freq, mean(bias(j).kc_ch4,2),'-');
   plot(calc.freq, mean(bias(j).sar_ch4,2),'-');
   %%%plot(calc.freq, mean(bias(j).sar_ch4 - bias.kc_ch4,2),'-');
   plot(calc.freq, stdv(j).kc_ch4,'-',calc.freq,stdv(j).sar_ch4,'-');
   %plot(calcs.freq, std(bias.kc_ch4 - bias.sar_ch4,0,2),'-')
   legend('kCARTA mean','SARTA mean','kC std.dev','SAR std.dev','location','best');
   ylabel('perturb minus unpert (K)');xlabel('wvn cm^{-1}');
   th1=title([csens ' Mean d(BT) from 10% CH4.pert']);set(th1,'Interpreter','none');
      %saveas(gcf, [phome 'kc_sar_pert_10pc_ch4_0deg.fig'],'fig');
      
% SO2
fign = fign + 1;
for j = 1:7
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.so2(j:8:392); % nadir
  bias(j).kc_so2  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_so2 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_so2  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_so2 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
 figure(fign);clf;hold on; grid on; xlim([1040 2600]);
   plot(kfreq, mean(bias(j).kc_so2,2),'-');
   plot(sfreq, mean(bias(j).sar_so2,2),'-');
   %plot(calc.freq, std(bias.kc_so2 - bias.sar_so2,0,2),'-') 
   plot(kfreq, stdv(j).kc_so2,'-',sfreq,stdv(j).sar_so2,'-');
   plot(sfreq(ix), mean(bias(j).sar_so2(ix,:),2) - mean(bias(j).kc_so2(iy,:),2),'-');
   legend('kCARTA mean','SARTA mean','kC std.dev','sar std.dev','location','best'); 
   xlabel('wvn cm^{-1}');ylabel('perturb minus unpert (K)');
   th1=title([csens ' Mean,std,bias d(BT) 10% so2.pert']);set(th1,'Interpreter','none');
   %saveas(gcf, [phome 'kc_sar_pert_10pc_so2_0deg.fig'],'fig');
   
% HNO3
fign = fign + 1;
for j = 1:7
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.hno3(j:8:392); % nadir
  bias(j).kc_hno3  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_hno3 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_hno3  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_hno3 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
 figure(fign);clf;grid on; hold on; xlim([700 1800]);
   plot(kfreq,   mean(bias(j).kc_hno3,2),'-');
   plot(sfreq,   mean(bias(j).sar_hno3,2),'-'); 
   plot(kfreq, stdv(j).kc_hno3,'-',sfreq,stdv(j).sar_hno3,'-');
   %plot(calc.freq(ix), std(bias.kc_hno3(ix) - bias.sar_hno3(ix),0,2),'-')
   plot(sfreq(ix), mean(bias(j).sar_so2(ix,:),2) - mean(bias(j).kc_so2(iy,:),2),'-');
   th1=title([csens ' Mean,std,bias d(BT) 10% hno3.pert']);set(th1,'Interpreter','none'); 
   legend('kCARTA mean','SARTA mean','kc std.dev','sar std.dev');
   xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
    %saveas(gcf, [phome 'kc_sar_pert_10pc_hno3_0deg.fig'],'fig');

% HN3
fign = fign + 1;
for j=1:7
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.nh3(j:8:392); 
  bias(j).kc_nh3  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_nh3 = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_nh3  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_nh3 = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
 figure(fign);clf;hold on; grid on; xlim([700 1200])
   plot(kfreq, mean(bias(j).kc_nh3,2),'-');
   plot(sfreq, mean(bias(j).sar_nh3,2),'-'); 
   plot(kfreq, stdv(j).kc_nh3,'-', sfreq,stdv(j).sar_nh3,'-');
   %plot(calc.freq(ix), mean(bias.kc_nh3(ix) - bias.sar_nh3(ix),2),'-')
   plot(sfreq(ix), mean(bias(j).kc_nh3(iy) - bias(j).sar_nh3(ix),2),'-');
   grid on; legend('kCARTA mean','SARTA mean','kCARTA std.dev','SARTA std.dev');
   xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
   th1=title([csens ' Mean,std,bias d(BT) 10% NH3.pert']);set(th1,'Interpreter','none'); 
   %saveas(gcf, [phome 'kc_sar_pert_10pc_nh3_0deg.fig'],'fig');

% CO
fign = fign + 1;
for j=1:7
  izo{j} = indx.unp(j:8:392);  xzo{j} = indx.co(j:8:392);
  bias(j).kc_co  = calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j});
  bias(j).sar_co = calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j});
  stdv(j).kc_co  = std(calc.bkc(:,xzo{j}) - calc.bkc(:,izo{j}),0,2);
  stdv(j).sar_co = std(calc.bsc(:,xzo{j}) - calc.bsc(:,izo{j}),0,2);
end
j=1;
 figure(fign);clf;plot(kfreq, mean(calc.bkc(:,cell2mat(izo)),2),'-')
      hold on; plot(sfreq, mean(calc.bsc(:,cell2mat(izo)),2),'-')
 figure(fign);clf;hold on; grid on; xlim([2000 2280])
   plot(kfreq, mean(bias(j).kc_co,2),'-');
   plot(sfreq, mean(bias(j).sar_co,2),'-'); 
   plot(kfreq, stdv(j).kc_co,'-',sfreq,stdv(j).sar_co,'-');
   %plot(sfreq, std(bias.kc_co - bias.sar_co,0,2),'-')
   plot(sfreq(ix), mean(bias(j).kc_co(iy) - bias(j).sar_co(ix),2),'-');
   title([csens ' Mean,std,bias d(BT) 10% CO.pert'],'Interpreter','none'); 
   legend('kCARTA mean','SARTA mean','kCARTA std','SARTA std','bias');
   xlabel('wvn cm^{-1}'); ylabel('perturb minus unpert (K)')
   %saveas(gcf, [phome 'kc_sar_pert_10pc_co_0deg.fig'],'fig');


% --------------------------------------------------------
% BIAS kc - sar unperturbed at the different secant angles 
% --------------------------------------------------------
fign = fign+1;

for j = 1:7
  izo{j} = indx.unp(j:8:392);
  bias(j).unp =  mean(calc.bkc(:,izo{j}) - calc.bsc(:,izo{j}), 2);
  stdv(j).unp =  std(calc.bkc(:,izo{j}) - calc.bsc(:,izo{j}), 0,2);
end

%bias.all(i,:) = mean(calcs.bkc(:,ja(i,:)) - calcs.bsc(:,ja(i,:)),2);

figure(fign);clf;hold on; 
  for j=1:2
    plot(calc.freq, bias(j).unp,'-');
  end

%}
