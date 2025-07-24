%
cd /home/chepplew/projects/sarta/prod_2019/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt, int2bits, mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/chepplew/projects/sarta/matlabcode

% Original R49 RTP:
srcdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
srcfn = [srcdr  'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];

[head hatt prof patt] = rtpread(srcfn);

% hardwire the kCARTA production run and coefficient calcs (matches existing 
%  directory paths)
prod_run = 'prod_2019';
buildver = 'dec2018';

% Choose which regression set to compare.
all_cregr = {'r49','saf704'};
cregr     = 'r49'; %'saf704';

% Enter which sensor to compare w/kcarta {'IASI','AIRS','CRIS'}
allsens = {'AIRS','CRIS_HR','IASI'};
csens = 'IASI';

% home for plots
phome = ['/home/chepplew/projects/sarta/' prod_run '/' lower(csens) '/figs/'];

% Set output directory and sarta calc rtp file name
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) ...
         '/' buildver '/tests/'];

% run SARTA w/ HDO

    x      = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq   = x.f_iasi;
    idchan = x.ichan_iasi;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
%    head.pmax  = 1013;        % check the kCARTA prediction surface pressure

    rtpx = [outdr 'r49_1013_400p_seaemis_nadir_8461.rtp'];
    tmp = mktemp();
    outfiles = rtpwrite_12(tmp,head,hatt,prof,patt);

    %SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_basic';
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_dec2018_hdo';

if(exist('/home/chepplew/logs/sar_out2.log'))
  delete('/home/chepplew/logs/sar_out2.log')
end

% -----------------------------------------------------
% run SARTA and read in results (save data if needed)
% -----------------------------------------------------

%switch csens
%  case 'IASI'

    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];

    eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sar_out2.log']);
    eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);

    cfin = [tmp '.sar'];
    [~,~,ptemp,~] = rtpread_12(cfin);

    % Save the calculations
    %rtpwrite_12(sarout,head,hatt,ptemp,patt);

    calc.bsc_nom  = rad2bt(head.vchan, ptemp.rcalc);

% --------------------------------------------------
% run sarta w/perturbed HDO (10% depleted)
  eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sar_out2.log']);
  eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);

  cfin = [tmp '.sar'];
  [~,~,ptemp,~] = rtpread_12(cfin);
  calc.bsc_d10 = rad2bt(head.vchan, ptemp.rcalc);

% run sarta w/perturbed HDO (10% enhanced)
  eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sar_out2.log']);
  eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);

  cfin = [tmp '.sar'];
  [~,~,ptemp,~] = rtpread_12(cfin);
  calc.bsc_e10 = rad2bt(head.vchan, ptemp.rcalc);

% --------------------------------------------
% Plot difference
% --------------------------------------------
figure(3);clf;
  plot(head.vchan, mean(calc.bsc_nom,2) - mean(calc.bsc_d10,2), '-'); 
    grid on; axis([1100 2800 -0.4 0.4]);title([csens ' SARTA.hdo Orig - 0.1 depleted'])
    xlabel('wavenumber cm^{-1}'); ylabel('dBT (K)');
    hold on; plot(head.vchan, std(calc.bsc_nom - calc.bsc_d10,0,2),'-');
    legend('mean','std')

% ================================================================
% compare w/ kCARTA
% /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/ ...
%    REGR49_400ppm_H2016_Dec2018_AIRS2834/hdo_radianceeffect.m
% ================================================================
wrkdir=pwd;
cd /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Dec2018_AIRS2834/
clear rad0* t0*
for ii = 1 : 49
  rad0   = load(['RAD1100_seaemis/convolved_kcarta_RAD1100_' num2str(ii) '_radiances.mat']);
  rad0p1 = load(['RAD1100_depleted0.1/convolved_kcarta_RAD1100_' num2str(ii) '_radiances.mat']);
  rad0p6 = load(['RAD1100_depleted0.6/convolved_kcarta_RAD1100_' num2str(ii) '_radiances.mat']);

  t0(ii,:,:)   = rad2bt(rad0.fiasi,rad0.riasi_all);
  t0p1(ii,:,:) = rad2bt(rad0.fiasi,rad0p1.riasi_all);
  t0p6(ii,:,:) = rad2bt(rad0.fiasi,rad0p6.riasi_all);
end
freq = rad0.fiasi;
nchan = length(freq);
cd(wrkdir)

t0   = permute(t0,[2 1 3]);
t0p1 = permute(t0p1,[2 1 3]);
t0p6 = permute(t0p6,[2 1 3]);

%plot(freq,rad2bt(freq,rad0.riasi_all));
%plot(freq,reshape(t0,nchan,49*8))

%plot(freq,reshape(t0,nchan,49*8)-reshape(t0p1,nchan,49*8))
figure(4);clf;
  plot(freq,mean((reshape(t0,nchan,49*8)-reshape(t0p1,nchan,49*8))'),'b',...
     freq,std((reshape(t0,nchan,49*8)-reshape(t0p1,nchan,49*8))'),'c--')
ylabel('\deltaBT'); xlabel('Wavenumber cm-1');
title('Orig - 0.1 depleted HDO');
hl = legend('mean','std','location','best');
grid on;

bias.kc_nad  = (reshape(t0(:,:,1),nchan,49)-reshape(t0p1(:,:,1),nchan,49));
bias.sar_nad = calc.bsc_nom - calc.bsc_d10;
stdv.kc_nad  = std(reshape(t0(:,:,1),nchan,49)-reshape(t0p1(:,:,1),nchan,49),0,2);
stdv.sar_nad = std(calc.bsc_nom - calc.bsc_d10,0,2);

figure(5);clf;
  plot(freq, mean(bias.kc_nad,2),'-', freq,mean(bias.sar_nad,2),'-' );
  grid on; ylabel('\deltaBT'); xlabel('Wavenumber cm-1'); 
  legend('kCARTA','SARTA');title('IASI kC, SAR Orig - 0.1 depl HDO');
  axis([1100 2800 -0.5 0.3])   
hold on;
  plot(freq, stdv.kc_nad,'-', freq, stdv.sar_nad,'-');
  plot(freq, std(bias.kc_nad - bias.sar_nad,0,2));
  legend('kCARTA mean','SARTA mean', 'kCARTA std','SARTA std','std kc - sar');
  %saveas(gcf, [phome 'kc_sar_iasi_hdo_pert_validation.fig'],'fig')
  
