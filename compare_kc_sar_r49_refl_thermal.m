% compare_kc_sar_r49_refl_thermal.m

cd /home/chepplew/projects/sarta/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt, int2bits, mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/chepplew/projects/sarta/matlabcode

% PreDefine configuration: prod_run, buildver, csens, cregr:
prod_run = 'prod_2021';
buildver = 'may2021';
csens    = 'CRIS_LR';
cregr    = 'R49'; %'saf704';

% Choose which regression set to compare.
all_cregr = {'R49','saf704'};

% Enter which sensor to compare w/kcarta {'IASI','AIRS','CRIS'}
allsens = {'AIRS_L1C','CRIS_LR','CRIS_HR','IASI','CHIRP'};

% home for plots
phome = ['/home/chepplew/projects/sarta/' prod_run '/' lower(csens) ...
         '/' buildver '/figs/'];

% Set output directory and sarta calc rtp file name
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) ...
         '/' buildver '/tests/'];

% -------------------------------------------------------------------------------
%      Set up RTP file for different surface properties (if not already)

outgn = ['/home/chepplew/data/sarta/' prod_run '/generic/'];
% fortp = 'r49_1013_400p_8angs_sfc_pert_v2.rtp';        % v2: incl. e 0.85 
% fortp = 'r49_1100_400p_8angs_sfc_pert_815zobs';
fortp  = 'r49_1013_400p_8angs_sfc_pert_815zobs.rtp';
fortp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert_v2.rtp';

if(~exist([outgn fortp]))    
  % Original R49 RTP:
  srcdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
  srcfn = [srcdr  'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];

  [hd ha pd pa] = rtpread(srcfn);

  [nr np] = size(pd.emis);
  % Different sensitivty tests: 
  % 1. decrease ocean emissivity by factor 4. to be more reflective
  emis.per2 = 1. - 4.*(1. - pd.emis(:,1));

  %2. Set a uniform 0.75 and 0.85 emissivity:
  emis.per3 = 0.85* ones(nr, 1);
  emis.per4 = 0.75* ones(nr, 1);

  %3. Set unity sfc emissivity
  emis.unit = ones(nr,1);

  % 4 extreme land from paper:
  emis.per5 = [0.93, 0.90, 0.91, 0.90, 0.90, 0.90, 0.85, ...
               0.79, 0.65, 0.70, 0.69, 0.78, 0.65, 0.58, ...
	       0.57, 0.58, 0.60, 0.62, 0.67]';
	     
  % Now duplicate the 49 profiles for different viewing angles:
  % SARTA limited to angles < 63-deg
  % *** ONLY do this once to set it up! ***
  angles = [0        8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
            65.0428 69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];

  % Do the nadir duplication and force rho = (1-e)/pi
  h2 = struct;
  p2 = struct;
  for ii = 1 : 49
     [hy,py] = replicate_rtp_headprof(hd,pd,ii,6);
     py.emis(:,1)   = pd.emis(:,1);
     py.emis(:,2)   = emis.per2;
     py.emis(:,3)   = emis.per3;
     py.emis(:,4)   = emis.per4;
     py.emis(:,5)   = emis.per5;
     py.emis(:,6)   = emis.unit;
     py.rho(:,2)    = (1.0 - py.emis(:,2))./3.14159;
     py.rho(:,3)    = (1.0 - py.emis(:,3))./3.14159;
     py.rho(:,4)    = (1.0 - py.emis(:,4))./3.14159;
     py.rho(:,5)    = (1.0 - py.emis(:,5))./3.14159;
     py.rho(:,6)    = (1.0 - py.emis(:,6))./3.14159;  
     if ii == 1
       h2 = hy;
       p2 = py;
     else
       [h2,p2] = cat_rtp(h2,p2,hy,py);
     end
     fprintf(1,'.')
  end
  clear hy py;
  fprintf(1,'\n')

  % Do the view angle duplications
  h3 = struct;
  p3 = struct;
  nprof = size(p2.satzen,2);     % expecting 294
  for ii = 1 : nprof
     [hy,py] = replicate_rtp_headprof(h2,p2,ii,7);  % was 8
     for jj = 2:7
       py.satzen(jj) = angles(jj);
     end
     if ii == 1
       h3 = hy;
       p3 = py;
     else
       [h3,p3] = cat_rtp(h3,p3,hy,py);
     end
     if(~mod(ii,49)); fprintf(1,'.'); end
  end
  clear hy py;
  fprintf(1,'\n')

  % write out pertubation set:
  rtpwrite([outgn fortp], h3, hatt, p3, patt);

end        % end: if test file does not exist
% ----------------------------------------------------------------------
% Prepare for the SARTA runs
%
% Algorithm for indexing each perturbation
% ----------------------------------------

indx = struct;
for ii=1:49
  for kk=1:8
    jj=8*(ii-1)+kk;
    %indx.emis1(jj)  = [1:49:392];
    indx.emis1(jj)  = kk  + (ii-1)*48;
    indx.emis2(jj)  = 8+kk + (ii-1)*48;
    indx.emis3(jj)  = 16+kk + (ii-1)*48;
    indx.emis4(jj)  = 24+kk + (ii-1)*48;
    indx.emis5(jj)  = 32+kk + (ii-1)*48;
    indx.emis6(jj)  = 40+kk + (ii-1)*48;
  end
end

% Reload RTP to do the test:
fnrtp = [outgn fortp];
[head hatt prof patt] = rtpread(fnrtp);
if(isfield(prof,'rcalc')) 
  prof = rmfield(prof,'rcalc');
end

% update header for chosen sensor and write rtp for SARTA 
% and select SARTA executable

switch csens
  case 'IASI'
    x      = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq   = x.f_iasi;
    idchan = x.ichan_iasi;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
      szz = size(prof.zobs);
      prof.zobs = 815000.0*ones(1,szz(2));
%    h2.pmax  = 1013;        % check the kCARTA prediction surface pressure
%    rtpx = [outdr 'r49_1013_400p_emis1_7angs_8461.rtp'];
    tmp = mktemp();
    outfiles = rtpwrite_12(tmp,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic';
    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic_optr';
    %SARTAEXE = '/asl/packages/sartaV108/Bin/sarta_iasi_may09_wcon_nte';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_thrm';

  case 'AIRS_L1C'
    hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf')
    freq  = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    %freq    = load('/home/chepplew/myLib/data/airs_f_l1c_lls.txt');
    %idchan  = int32([1:2645])';
    head.vchan = single(freq);
    head.ichan = int32(idchan);
    head.nchan = length(idchan);
      szz = size(prof.zobs);
      prof.zobs = 705000.0*ones(1,szz(2));
    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);

    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19';
    %%SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_2834_mar19_basic_optr_tra_nte_thrm';
    %%SARTAEXE = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';  

  case 'CRIS_LR'
    x       = load('/home/chepplew/myLib/data/fcris_lowres_4grd.mat');
    freq    = x.vchan;
    idchan  = x.idchan;
    head.vchan = freq';
    head.ichan = idchan';
    head.nchan = length(idchan);
    %head.pmax  = 1013;    
    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_lrg4_p2021_dev';

  case 'CRIS_HR'
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;    
    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_basic_optr_co2';

  case 'CHIRP'
    x = load('/home/chepplew/myLib/data/chirp_1702_freq.mat');
    head.vchan = single(x.vchan);
    head.ichan = int32([1:length(x.vchan)]');
    head.nchan = length(x.vchan);
    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_base_tra_thrm';

end
junk   = strsplit(SARTAEXE,'/');
fnout  = ['sar_' junk{end} '_' fortp];
sarout = ['outd/' fnout];

% --------------------------------
% Get the SARTA predictions
% --------------------------------
if(exist('/home/chepplew/logs/sarta/sar_out.log'))
  delete('/home/chepplew/logs/sarta/sar_out.log')
end

switch csens
  case 'IASI'
    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];
    eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sarta/sar_out1.log']);
    eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);
    cfin = [tmp '.sar'];
    [~,~,ptemp,~] = rtpread_12(cfin);
    calc.bsc  = rad2bt(head.vchan, ptemp.rcalc);

    % Save the calculations
    %rtpwrite_12(sarout,head,hatt,ptemp,patt);

  case 'AIRS_L1C'
    ifn = tmp;      ofn = [tmp '.sar'];
    eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    %command=[SARTAEXE ' fin=rtpx fout=' sarout ' > /home/chepplew/logs/sarta/sar_out.log'];
    %system(command);
    [~,~,ptemp,~] = rtpread(ofn);
    calc.bsc = rad2bt(head.vchan, ptemp.rcalc);

  case 'CRIS_LR'
    ifn = tmp;      ofn = [tmp '.sar'];
    tic;
      eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    toc
    %eval(['! ' SARTAEXE ' fin=rtpx  fout=' sarout]);
    [~,~,ptemp,~] = rtpread(ofn);
    calc.rad = ptemp.rcalc;
    calc.bsc = rad2bt(head.vchan, ptemp.rcalc);

  case 'CRIS_HR'
    ifn = tmp;      ofn = [tmp '.sar'];
    tic;
      eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
    toc
    %eval(['! ' SARTAEXE ' fin=rtpx  fout=' sarout]);
    [~,~,ptemp,~] = rtpread(ofn);
    calc.bsc = rad2bt(head.vchan, ptemp.rcalc);

  case 'CHIRP'
    ifn = tmp;      ofn = [tmp '.sar'];
    command=[SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log'];
    tic
      system(command);
    toc
    [~, ~, ptemp, ~] = rtpread(ofn);
    calc.bsc = rad2bt(head.vchan, ptemp.rcalc);

    % Save the calculations
end

% ------------------------------------------
% Get the kCARTA TOA rads
% ------------------------------------------
switch prod_run
  case 'prod_2019'
    %kpath = '/home/chepplew/data/kcarta/REGR49/400p_1013mb_emis_p75_7angs/';
    %kpath = '/home/chepplew/data/kcarta/REGR49/400p_1013mb_emis1_7angs/';
    khome = '/home/chepplew/data/sarta/prod_2019/generic/kcarta/';
    %kpath = [khome 'REGR49/surface_pert/'];
    %kpath = [khome 'REGR49/surface_pert_v2/'];
    %kpath = [khome 'REGR49/surface_pert_CKD25/'];
    kpath = [khome 'REGR49/r49_1013_400p_8angs_sfc_pert_815zobs/'];
    klist = dir([kpath 'individual_prof_convolved_kcarta_crisHI_*.mat']);
  case 'prod_2020' 
    kpath = ['/home/chepplew/data/sarta/prod_2020/generic/kcarta/' ...
             '400p_8angs_sfc_pert_815zobs/'];
    % r49_1013_400p_8angs_sfc_pert_815zobs
    klist = dir([kpath 'convolved_kcarta_RADTEST_SARTA_ChrisH_*.mat']);
  case 'prod_2021'
    kpath = ['/home/chepplew/data/sarta/prod_2021/generic/kcarta/REGR49/' ...
             '400p_1013mb_8angs_sfc_pert_v2/'];
    klist = dir([kpath 'individual_prof_convolved_kcarta*.mat']);  

end
whos klist

clear prfnums fnparts;
% reorder these in profile number order
for i=1:length(klist) 
  fnparts{i} = strsplit(klist(i).name, {'_' '.'}); 
  prfnums(i) = str2num(cell2mat(fnparts{i}(7)));       % 6 or 7
end
[IB IC] = sort(prfnums);  
krc = []; 
for i=IC 
  x   = load([klist(i).folder '/' klist(i).name]);
  switch csens
    case 'IASI'
      krc    = [krc x.riasi_all]; 
      if(i == 1) kc_frq = x.fiasi; end
    case 'CRIS_LR'
      krc    = [krc x.lo_rcris_all];
      if(i == 1) kc_frq = x.lo_fcris; end
    case 'CRIS_HR'
      krc    = [krc x.hi_rcris_all];
      if(i == 1) kc_frq = x.hi_fcris; end
    case 'AIRS_L1C'
      krc    = [krc x.rKc];
      if(i == 1) kc_frq = x.fKc; end
    case 'CHIRP'
      krc   = [krc x.med_rcris_all];
      if(i == 1) kc_frq = x.med_fcris; end
  end
end

calc.rkc = krc;
calc.bkc = rad2bt(kc_frq, krc);
calc.frq = head.vchan;
calc.rks = krc - ptemp.rcalc;
clear krc;

% Stats for sea emissivity
radstd   = nanstd(ptemp.rcalc(:,indx.emis1) - krc(:,indx.emis1),0,2);
 cdbm    = 0.5*( nanmean(calc.bsc(:,indx.emis1),2) + nanmean(calc.bkc(:,indx.emis1),2) );
 mdr     = 1E-3*( 1./drdbt(calc.frq,cdbm) );
btstd    = mdr.*radstd;
btser    = btstd./sqrt(size(indx.emis1,2));

%{
% -----------------------------------------------
%              PLOTTING SECTION
% -----------------------------------------------
if(strcmp(csens,'AIRS_L1C')) 
  [~, ib] = sort(head.vchan);
else
  ib = find(head.vchan);
end

% Display SFC emissivity used
fign = 1;
figure(fign); clf; hold on;
  plot(ptemp.efreq(:,indx.emis1(1)), ptemp.emis(:,indx.emis1(1)),'.-')
  plot(ptemp.efreq(:,indx.emis2(1)), ptemp.emis(:,indx.emis2(1)),'.-')
  plot(ptemp.efreq(:,indx.emis3(1)), ptemp.emis(:,indx.emis3(1)),'.-')
  plot(ptemp.efreq(:,indx.emis4(1)), ptemp.emis(:,indx.emis4(1)),'.-')
  plot(ptemp.efreq(:,indx.emis5(1)), ptemp.emis(:,indx.emis5(1)),'.-')
  plot(ptemp.efreq(:,indx.emis6(1)), ptemp.emis(:,indx.emis6(1)),'.-')
  legend('emis1','emis2','emis3','emis4','emis5','emis6');
  grid on;ylim([0.45 1.05]);xlabel('wvn cm^{-1}')
  title('surface emissivity variants')
  %saveas(gcf, [phome 'sfc_emissivity_variants.png'],'png');

% Hack to get all profiles for emisX for first 6 angles
iwnt = [];
for j = 1:6
  iwnt = [iwnt indx.emis1(j:8:392)];
end
iwnt = sort(iwnt);

% Stats of Selected angles & emissivity of SARTA - kCARTA
bias = mean(calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt),2);
stdv = std(calc.bsc(ib,iwnt) - calc.bkc(ib,iwnt),0,2);
freq = head.vchan(ib);
fign=2;
figure(fign);clf; hold on; grid on;
  plot(freq,bias,'-', freq, stdv,'-');
  plot(head.vchan(ib), mean(calc.bsc(ib,indx.emis1(iwnt)),2),'b-');
  plot(head.vchan(ib), mean(calc.bsc(ib,indx.emis2(iwnt)),2),'c-'); 
  legend('sar.e1','sar.e2')
  
fign = fign+1;
figure(fign);clf;hold on; for j=1:7
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis1(j:8:392)),2) - ...
       mean(calc.bsc(ib,indx.emis1(j:8:392)),2),'-'); end
  title('kCARTA')

fign = fign+1;
th3 = 'kC vs SARTA'; 
fh3=figure(fign);clf; %set(fh3,'resize','off');set(fh3,'Position',fh3.Position+[0 0 0 120]);
 h1=subplot(211);hold on; j=1; grid on;
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis1(iwnt)),2),'b-');
  plot(head.vchan(ib), mean(calc.bsc(ib,indx.emis1(iwnt)),2),'c-');
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis4(iwnt)),2),'r-');
  plot(head.vchan(ib), mean(calc.bsc(ib,indx.emis4(iwnt)),2),'m-');
 h2=subplot(212);hold on;grid on;
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis1(iwnt)),2) - ...
    mean(calc.bsc(ib,indx.emis1(iwnt)),2),'-');
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis3(iwnt)),2) - ...
    mean(calc.bsc(ib,indx.emis3(iwnt)),2),'-');
  plot(head.vchan(ib), mean(calc.bkc(ib,indx.emis4(iwnt)),2) - ...
    mean(calc.bsc(ib,indx.emis4(iwnt)),2),'-');
  legend('sea','0.85','0.75'); ylim([-0.4 0.4]);  
  title('SARTA')

% ============================================================
% Reflected component relative to e=1 separately for SARTA and kCARTA 
for j=1:7
  bias.sar_1(:,j) = mean(calc.bsc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bsc(:,indx.emis1(j:8:392)),2);
  bias.kc_1(:,j)  = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bkc(:,indx.emis1(j:8:392)),2);
  bias.sar_2(:,j) = mean(calc.bsc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bsc(:,indx.emis2(j:8:392)),2);
  bias.kc_2(:,j)  = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bkc(:,indx.emis2(j:8:392)),2);
  bias.sar_3(:,j) = mean(calc.bsc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bsc(:,indx.emis3(j:8:392)),2);
  bias.kc_3(:,j)  = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bkc(:,indx.emis3(j:8:392)),2);
  bias.sar_4(:,j) = mean(calc.bsc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bsc(:,indx.emis4(j:8:392)),2);
  bias.kc_4(:,j)  = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bkc(:,indx.emis4(j:8:392)),2);
  bias.sar_5(:,j) = mean(calc.bsc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bsc(:,indx.emis5(j:8:392)),2);
  bias.kc_5(:,j)  = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                  - mean(calc.bkc(:,indx.emis5(j:8:392)),2);
end

j = 1; % (j=1: 0-deg. j=7: 60-deg)
fign = fign + 1;
figure(fign); clf(fign); hold on;
  plot(calc.frq, bias.sar_1(:,j),'-', calc.frq, bias.kc_1(:,j),'-');
  plot(calc.frq, bias.sar_2(:,j),'-', calc.frq, bias.kc_2(:,j),'-');
  plot(calc.frq, bias.sar_3(:,j),'-', calc.frq, bias.kc_3(:,j),'-');
 grid on;xlim([600 2600])
 xlabel('wvn cm-1');ylabel('dBT (K)');
 legend('sea.sar','sea.kc','0.25sea.sar','0.25sea.kc','0.85.sar','0.85.kc')  
 title('Refl Thermal SARTA, kCARTA rel. to e=1, 0-deg')

% Bias between SARTA and kCARTA. All angles
% emis1: sea, emis2: sea*0.25, emis3: 0.85, emis4: 0.75, emis5: land, emis6:1
for j=1:7
  calc.bias_e1(:,j) = mean(calc.bkc(:,indx.emis1(j:8:392)),2) ...
                   - mean(calc.bsc(:,indx.emis1(j:8:392)),2);
  calc.bias_e3(:,j) = mean(calc.bkc(:,indx.emis3(j:8:392)),2) ...
                   - mean(calc.bsc(:,indx.emis3(j:8:392)),2);
  calc.bias_e4(:,j) = mean(calc.bkc(:,indx.emis4(j:8:392)),2) ...
                   - mean(calc.bsc(:,indx.emis4(j:8:392)),2);
  calc.bias_e6(:,j) = mean(calc.bkc(:,indx.emis6(j:8:392)),2) ...
                   - mean(calc.bsc(:,indx.emis6(j:8:392)),2);
  calc.stdv_e1(:,j) = std(calc.bkc(:,indx.emis1(j:8:392)) ...
                   - calc.bsc(:,indx.emis1(j:8:392)),0,2);
  calc.stdv_e3(:,j) = std(calc.bkc(:,indx.emis3(j:8:392)) ...
                   - calc.bsc(:,indx.emis3(j:8:392)),0,2);
  calc.stdv_e4(:,j) = std(calc.bkc(:,indx.emis4(j:8:392)) ...
                   - calc.bsc(:,indx.emis4(j:8:392)),0,2);
  calc.stdv_e6(:,j) = std(calc.bkc(:,indx.emis6(j:8:392)) ...
                   - calc.bsc(:,indx.emis6(j:8:392)),0,2);
end

fign = fign+1;
fh4=figure(fign);clf; set(fh4,'resize','off');set(fh4,'Position',fh4.Position+[0 0 0 180]);
  h41=subplot(411); hold on;grid on;
  for j=1:6 
    plot(head.vchan(ib), calc.bias_e1(ib,j),'-'); 
    plot(head.vchan(ib), calc.stdv_e1(ib,j),':'); 
  end
  ylim([-1.2 1.4]);title('kC - SAR sea emissivity vs satzen');ylabel('dBT (K)')  
  h42=subplot(412); hold on;grid on;
  for j=1:6 
    plot(head.vchan(ib), calc.bias_e3(ib,j),'-'); 
    plot(head.vchan(ib), calc.stdv_e3(ib,j),':'); 
  end
  ylim([-1.2 1.4]);title('kc - SAR 0.85 emissivity')  
  h43=subplot(413); hold on;grid on;
  for j=1:6 
    plot(head.vchan(ib), calc.bias_e6(ib,j),'-'); 
    plot(head.vchan(ib), calc.stdv_e6(ib,j),':'); 
  end
  ylim([-1.2 1.4]);title('kC - SAR unit emissivity');xlabel('wavenumber cm^{-1}'); 
  h44=subplot(414); hold on;grid on;
  for j=1:6 
    plot(head.vchan(ib), calc.bias_e4(ib,j),'-'); 
    plot(head.vchan(ib), calc.stdv_e4(ib,j),':'); 
  end
  ylim([-1.2 1.4]);title('kC - SAR 0.75 emissivity');xlabel('wavenumber cm^{-1}'); 
  %saveas(gcf,[phome 'cris_lr_kc_sar_r49_reflthrm_mean_stdv_vs_emiss_satzen.fig'],'fig');
    
%}
