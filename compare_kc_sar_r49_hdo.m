% compare_kc_sar_r49_hdo.m
%
% Compute simulations on 49 regression profiles
% with different HDO depletions using  prof.udef(20,:)
%   Depletion in pp.mil. (e.g. typically 0 to 800)

addpath /asl/matlib/h4tools                      % rtpread
addpath /asl/matlib/aslutil                      %  rad2bt

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


% Original Test Profiles
if(contains(opts.csens,'PBL'))
  srcrtp = ['/home/chepplew/data/sarta/prod_2025/generic/regr49_pbl_nh3_op.rtp'];
else
   srcrtp = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/' ...
             'RUN_KCARTA/REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];
   srcrtp = '/home/chepplew/gitLib/sarta/test/regr49_1013_400ppm.op.rtp';
end

% --------------------------
% kCARTA source
% --------------------------
kpath = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
        'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];

klist.hdo = dir([kpath 'RAD1013_depleted0.6/convolved_kcarta_RAD1013_*_radiances.mat']);
klist.ref = dir([kpath 'RAD1013/convolved_kcarta_RAD1013_*_radiances.mat']);

whos klist

% sort kCARTA files in profile order:
clear prfnums fnparts;
% reorder these in profile number order
for i=1:length(klist.hdo)
  fnparts{i} = strsplit(klist.hdo(i).name, {'_' '.'});
  prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
end
[B IB] = sort(prfnums);
krc = [];
for i=IB
  x   = load([klist.hdo(i).folder '/' klist.hdo(i).name]);
  if(strcmp(opts.csens,'IASI'))     krc = [krc x.rKcIasi]; end
  if(strcmp(opts.csens,'CRIS_HR'))  krc = [krc x.rcris_all]; end
  if(strcmp(opts.csens,'AIRS_L1C')) krc = [krc x.rKc]; end
  if(strcmp(opts.csens,'CHIRP'))    krc = [krc x.rKc]; end
end
switch opts.csens
  case 'IASI'
    kc_frq = x.fiasi;
  case 'CRIS_HR'
    kc_frq = x.hi_fcris;
  case 'AIRS_L1C'
    kc_frq = x.fKc;
  case 'CHIRP'
    kc_frq = x.med_fcris;
end

rkc.hdo = krc;
bkc.hdo = rad2bt(kc_frq, krc);
freq    = kc_frq;

% Get  reference kCARTA TOA Rads
krc = [];
for i=IB
  x   = load([klist.ref(i).folder '/' klist.ref(i).name]);
  if(strcmp(csens,'IASI'))     krc = [krc x.rKcIasi]; end
  if(strcmp(csens,'CRIS_HR'))  krc = [krc x.rcris_all]; end
  if(strcmp(csens,'AIRS_L1C')) krc = [krc x.rKc]; end
end
rkc.ref = krc;
bkc.ref = rad2bt(kc_frq, krc);
clear krc;

%{
figure(2);clf;plot(freq, nanmean(bkc.hdo(:,[1:8:392]),2),'-')
  hold on; plot(freq, nanmean(bkc.ref(:,[1:8:392]),2),'-')
%}
% ================================================
%   GET SARTA rcalc
% ===============================================
% create temporary workspace
tempDir = mktemp();

% Load original test profiles
[head, hatt, prof, patt] = rtpread(srcrtp);

% set spreas = std. atm.
prof.spres = 1013.25*ones(size(prof.spres));

% add HDO depletion
prof.udef(20,:) = zeros(size(prof.rlat));

switch opts.csens
  case {'AIRS_L1C','AIRS_PBL'}
    fn.r49_op = [tempDir '_r49_' lower(csens) '_op.rtp'];
    load('/home/chepplew/myLib/data/airs_f.mat','fairs');
    load('/home/sbuczko1/git/rtp_prod2/airs/static/sarta_chans_for_l1c.mat','ichan');
    hinfo  = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf')
    freq   = hdfread(hinfo.SDS(2));
    idchan = hdfread(hinfo.SDS(1));
    head.vchan = single(freq(:));
    head.ichan = int32(idchan(:));
    head.nchan = length(head.ichan);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_v3';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_pbl_dev_debug'

  case {'IASI','IASI_PBL'}
  
  case {'CRIS_HR','CRIS_HR_PBL'}
    fn.r49_op = [tempDir '_r49_' lower(csens) '_op.rtp'];
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;    
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_may18_hdo_x';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/cris_hr_pbl_dev';

  case {'CHIRP','CHIRP_PBL'}
    fn.r49_op = [tempDir '_r49_chirp_pbl_op.rtp'];
    load('/home/chepplew/myLib/data/chirp_1702_freq_4grd.mat','vchan');
    load('/home/chepplew/myLib/data/chirp_1702_freq_4grd.mat','ichan');
%    x       = load('/home/chepplew/myLib/data/chirp_1691_g2.mat');
%    x       = load('/home/chepplew/myLib/data/chirp_1691_freq.mat');
    head.vchan = single(vchan(:));
    head.ichan = int32(ichan(:));
    head.nchan = length(vchan);
    rtpwrite(fn.r49_op, head,hatt,prof,patt);
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_feb20_base_so2_nh3_hno3_n2o';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_p2022jul22_xnte_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta_f90/bin/chirp_dev';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_reg_prod';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/chirp_pbl_p2025jan25a_dev';
 
end
if(strcmp(opts.csens,'IASI'))
    outfns = rtpwrite_12(fn.r49_op, head,hatt,prof,patt);
else
    rtpwrite(fn.r49_op, head,hatt,prof,patt);
end

% -------------------------
% Run SARTA
% ------------------------

disp('SARTA running')
switch csens
  case 'IASI'
  
  case {'CRIS_HR','AIRS_L1C','CHIRP','AIRS_PBL','IASI_PBL','CRIS_PBL','CHIRP_PBL'}
  
    ifn = fn.r49_op;      ofn = [tempDir '_sar_r49_' lower(csens) '.rtp'];
    command=[SARTAEXE ' fin=' ifn ' fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']; 

end

% Run again using depleted HDO
tic; system(command); toc 
[hds,~,ptemp_d00,~] = rtpread(ofn);
calc.rad_d00 = ptemp_d00.rcalc;
calc.bsc_d00 = real(rad2bt(hds.vchan, ptemp_d00.rcalc));
 
prof.udef(20,:) = 600.0*ones(size(prof.rlat));
rtpwrite(fn.r49_op,head,hatt,prof,patt)
tic; system(command); toc
[~,~,ptemp_d06,~] = rtpread(ofn);
calc.rad_d06 = ptemp_d06.rcalc;
calc.bsc_d06 = real(rad2bt(head.vchan, ptemp_d06.rcalc));

prof.udef(20,:) = 1200.0*ones(size(prof.rlat));
rtpwrite(fn.r49_op,head,hatt,prof,patt)
tic; system(command); toc
[~,~,ptemp_d12,~] = rtpread(ofn);
calc.rad_d12 = ptemp_d12.rcalc;
calc.bsc_d12 = rad2bt(head.vchan, ptemp_d12.rcalc);

prof.udef(20,:) = -600.0*ones(size(prof.rlat));
rtpwrite(fn.r49_op,head,hatt,prof,patt)
tic; system(command); toc
[~,~,ptemp_e06,~] = rtpread(ofn);
calc.rad_e06 = ptemp_e06.rcalc;
calc.bsc_e06 = real(rad2bt(head.vchan, ptemp_e06.rcalc));


%{
[sfreq, iss] = sort(head.vchan);
figure(2);clf;plot(freq, nanmean(calc.bsc_d06,2),'-')
  hold on; plot(freq, nanmean(calc.bsc_d6t2),'-')
clf; hold on; grid on;
  plot(sfreq, mean(calc.bsc_d06,2) - mean(calc.bsc_d00,2),'-')
  plot(sfreq, mean(calc.bsc_d12,2) - mean(calc.bsc_d00,2),'-')

%}




