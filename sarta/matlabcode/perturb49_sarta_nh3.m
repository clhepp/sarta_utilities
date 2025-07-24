% perturb49_sarta_nh3.m
%
% FOr use with CrIS HiRes SARTA with NH3
%     Take original 49 regression profiles and perturb NH3
% Compare with kCARTA predicts:
% /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/CO2_CO_coljac
% and NH3_coljac, SO2_coljac copied to: /home/chepplew/data/kcarta/sarta_breakouts/
% 
% Calulate TOA one pertubed gas at a time.
%
% NB: run by hand - selected sections.
%

cd /home/chepplew/projects/sarta/prod_2018/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                       % rad2bt

% Choose which SARTA to use (CrIS hiRes or AIRS L1b)
allvers = {'cris_hr','airs'};
myvers  = lower(myvers);
CRIS_HR=false;  AIRS=false;
if(~ismember(myvers,allvers)) error('invalid model version'); end
if(strcmp(myvers,allvers{1}))
  CRIS_HR  = true;
  csarta   = 'crisg4_may18_nh3';
  SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_may18_nh3'; end
if(strcmp(myvers,allvers{2}))
  AIRS     = true;
  csarta   = 'sarta_apr08_m140_wcon_nte_nh3';
  SARTAEXE = ['/home/chepplew/projects/sarta/sartaV108_juying/'...
              'NH3_sartaV108/BinV201/sarta_apr08_m140_wcon_nte_nh3'];
end  

% Minor gas list to perturb:
pert_cgas  = {'co','so2','nh3'};
pert_ngas  = [5, 9, 11];

% Original regression profile RTP file:
rtp49 = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP',...
         '/stdNH3_1100mb_op_400ppm.rtp'];
%{
% convert this to end at 1013 mb
[hd ha pd pa] = rtpread(rtp49);
varnames = fieldnames(pd);
varwant  = {'plevs','palts','ptemp','gas_1','gas_2','gas_3','gas_4','gas_5',...
            'gas_6','gas_9','gas_11','gas_12'};
zlevs=[98:101];
for k=1:numel(varwant)
    pd.(varwant{k})(zlevs,:) = 0.0;
end
pd.nlevs = 98*ones(size(pd.nlevs)); 
pd.spres = 1013.2*ones(size(pd.spres));
%}
	 
%rtp49 = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP',...
%         '/regr_rtp_6angs_49profs_1013mb_unitemis.rtp'];
[hd ha pd pa] = rtpread(rtp49);

% over-ride sfc emissivity = 1 or 0;
pd.emis = ones(size(pd.emis));

% Make header for CrIS hiRes w/ 4 guard channels (load hdr2) or AIRS L1b.
if CRIS_HR
  load('/home/chepplew/projects/sarta/prod_2016/cris_hr/rtp_head_str_2235g4.mat');
  hdr2.pfields = 1;
  hdr2.glist = sort([hdr2.glist; 11]);
  hdr2.gunit = [hdr2.gunit; 1];
  rtp49x = ['./sarta_data/' myvers '/r49_std_nh3_op_400_2235g4.rtp'];
end
if AIRS
  load('/home/chepplew/projects/airs/rtp_head_stru_airibrad_2378.mat','head')
  head.glist = sort([head.glist; 2; 4; 5; 6; 9; 11; 12]);
  head.gunit = ones(size(head.glist)); 
  head.ngas    = size(head.glist,1);
  head.ptype   = 1;
  head.pfields = 1;
  head.pmin    = 0.0050;
  hdr2   = head; 
  rtp49x = ['./sarta_data/' myvers '/r49_std_nh3_op_400_2378.rtp'];
end

if exist(rtp49x) disp(['Warning: ' rtp49x ' already exists.']); 
 ans = input('Overwrite (Y/N)?','s');
 if(strcmp(ans,'Y'))  rtpwrite(rtp49x,hdr2,ha,pd,pa); end
 if(strcmp(ans,'N')); end
else
 rtpwrite(rtp49x,hdr2,ha,pd,pa);
end

% duplicate original profiles amplifying selected gas (alone)

switch pert_ngas
  case 5
    k = 1;
    py = pd;
    py.gas_5  = 1.1* pd.gas_5;                % CO
    py.gas_2  = 1.1* pd.gas_2;                % CO2
    fortp = ['../sarta_data/' myvers '/r49_pert_' pert_cgas{k} '_400_' myvers '_op.rtp'];
    rtpwrite(fortp,hdr2,ha,py,pa);
    % Run sarta for the perturbed gas
    pfo = ['../sarta_data/' myvers '/sar_' csarta '_r49_pert_' pert_cgas{k} ...
           '_v03.rtp'];
    if exist('sym_pfo') delete('sym_pfo'); end
    unix([' ln -s ' pfo  ' sym_pfo '])
    command = [ SARTAEXE ' fin=' fortp ' fout=sym_pfo > ' ...
                  '/home/chepplew/logs/sarta/sar_stdout_pert.txt'];
    unix(command)

  case 9
    k = 2;
    py = pd;
    py.gas_9  = 1.1* pd.gas_9;                % SO2
    fortp = ['./sarta_data/' myvers '/r49_pert_' pert_cgas{k} '_400_' myvers '_op.rtp'];
    rtpwrite(fortp,hdr2,ha,py,pa);
    % Run sarta for the perturbed gas
    pfo = ['./sarta_data/' myvers '/sar_' csarta '_r49_pert_'  pert_cgas{k} ...
           '_v11.rtp'];
    if exist('sym_pfo') delete('sym_pfo'); end
    unix([' ln -s ' pfo  ' sym_pfo '])
    command = [ SARTAEXE ' fin=' fortp ' fout=sym_pfo > ' ...
                  '/home/chepplew/logs/sarta/sar_stdout_pert.txt'];
    unix(command)

  case 11
    k = 3;
    py = pd;
    py.gas_11 = 1.1* pd.gas_11;               % NH3
    fortp = ['./sarta_data/' myvers '/r49_pert_' pert_cgas{k} '_400_' myvers '_op.rtp'];
    rtpwrite(fortp,hdr2,ha,py,pa);
    % Run sarta for the perturbed gas
    pfo = ['./sarta_data/' myvers '/sar_' csarta '_r49_pert_100x_' pert_cgas{k} ...
           '_1100mb_unitemis_v01.rtp'];
    if exist('sym_pfo') delete('sym_pfo'); end
    unix([' ln -s ' pfo  ' sym_pfo '])
    command = [ SARTAEXE ' fin=' fortp ' fout=sym_pfo  > ' ...
                  '/home/chepplew/logs/sarta/sar_stdout_pert.txt'];
    unix(command)
end


% run sarta for the un-perturbed gas
fin = rtp49x;
ufo = ['./sarta_data/' myvers '/sar_' csarta '_r49_unpert_100xnh3_unitemis_v01.rtp'];
if exist('sym_ufo') delete('sym_ufo'); end
unix([' ln -s ' ufo  ' sym_ufo '])
command = [ SARTAEXE ' fin=' rtp49x ' fout=sym_ufo  > ' ...
                  '/home/chepplew/logs/sarta/sar_stdout_unp.txt'];
unix(command)

% Compare the two calculations
[hdu hau pdu pau] = rtpread('sym_ufo');
[hdp hap pdp pap] = rtpread('sym_pfo');

srtu  = rad2bt(hdu.vchan, pdu.rcalc);
srtum = nanmean(srtu,2);
srtp  = rad2bt(hdp.vchan, pdp.rcalc);
srtpm = nanmean(srtp,2);
srbias_mn = nanmean(srtp - srtu,2);
srbias_sd = nanstd(srtp - srtu,0,2);

% ------------------------------------------------------------
% compare w/ kcarta predicts (choose manually)
% ------------------------------------------------------------

nh3.kc_home = ['/home/chepplew/data/kcarta/sarta_breakouts/prod2018/',...
               'REGR49_400ppm_H2016_Mar2018_NH3/PERT/'];
nh3.kc_home = ['/home/chepplew/data/kcarta/sarta_breakouts/prod2018/'...
               'REGR49_400ppm_H2016_Mar2018_NH3/PERT/NH3_coljac_unitemiss/'];
nh3.kc_home = ['/home/chepplew/data/kcarta/sarta_breakouts/prod2018/'...
               'REGR49_400ppm_H2016_Mar2018_NH3/PERT/NH3_coljac_zeroemiss/'];

so2.kc_home = '/home/chepplew/data/kcarta/sarta_breakouts/SO2_coljac/';
co2.kc_home = '/home/chepplew/data/kcarta/sarta_breakouts/CO2_CO_coljac/';

% Load up kCARTA CO
% ------------------	   
co.kctu = zeros(49,2235);
co.kctp = zeros(49,2235);
% coljac col:1 10% CO jac. Col:2 +1 K change in entire T(z),  Col:3 stemp changed by 1 K
for i=1:49
  x = load([co2.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '.mat']);
  co.kctu(i,:) = rad2bt(x.fcris,x.rcris_all);

  y = load([co2.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '_coljac.mat']);
  co.kctp(i,:) = rad2bt(y.fcris,y.rcris_all(:,1));
  fprintf(1,'.');
end

co.kctum = nanmean(co.kctu,1);
co.kctpm = nanmean(co.kctp,1);
co.kctus = nanstd(co.kctu,0,1);
co.kctps = nanstd(co.kctp,0,1);
co.kcbias_mn = nanmean(co.kctp - co.kctu,1);
co.kcbias_sd = nanstd(co.kctp - co.kctu,0,1);

% Load up kCARTA NH3
% ------------------	   
nh3.kctu = zeros(49,2235);
nh3.kctp = zeros(49,2235);
% coljac col:1 10% NH3 jac. Col:2 +1 K change in entire T(z),  Col:3 stemp changed by 1 K
for i=1:49
  x = load([nh3.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '.mat']);
  nh3.kctu(i,:) = rad2bt(x.fcris,x.rcris_all);

  y = load([nh3.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '_coljac.mat']);
  nh3.kctp(i,:) = rad2bt(y.fcris,y.rcris_all(:,1));
  fprintf(1,'.');
end

nh3.fc    = x.fcris;
nh3.kctum = nanmean(nh3.kctu,1);
nh3.kctpm = nanmean(nh3.kctp,1);
nh3.kctus = nanstd(nh3.kctu,0,1);
nh3.kctps = nanstd(nh3.kctp,0,1);
nh3.kcbias_mn = nanmean(nh3.kctp - nh3.kctu,1);
nh3.kcbias_sd = nanstd(nh3.kctp - nh3.kctu,0,1);

% Load up kCARTA SO2
% ------------------	   
so2.kctu = zeros(49,2235);
so2.kctp = zeros(49,2235);
% coljac col:1 10% SO2 jac. Col:2 +1 K change in entire T(z),  Col:3 stemp changed by 1 K
for i=1:49
  x = load([so2.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '.mat']);
  so2.kctu(i,:) = rad2bt(x.fcris,x.rcris_all);

  y = load([so2.kc_home 'individual_prof_convolved_kcarta_AIRS_crisHI_' num2str(i) '_coljac.mat']);
  so2.kctp(i,:) = rad2bt(y.fcris,y.rcris_all(:,1));
  fprintf(1,'.');
end

so2.kctum = nanmean(so2.kctu,1);
so2.kctpm = nanmean(so2.kctp,1);
so2.kctus = nanstd(so2.kctu,0,1);
so2.kctps = nanstd(so2.kctp,0,1);
so2.kcbias_mn = nanmean(so2.kctp - so2.kctu,1);
so2.kcbias_sd = nanstd(so2.kctp - so2.kctu,0,1);

%{
frq = hdr2.vchan;

figno = 1;
figure(figno);plot(x.fcris,kctpm - kctum ,'-')     

nf9=figure(9);clf;plot(frq, srtpm,'-', frq, srtum,'-');grid on;
figure(9);clf;plot(frq, srtpm - srtum, '-'); grid on; xlim([650 1800]);
  xlabel('wavenumber cm-1');ylabel('pert minus unpert (K)');
  title('');
  
figure(2);clf;plot(fc, srtu(:,1) - srtp(:,1),'-');grid on;xlim([640 1500]);

nf10=figure(10);clf;plot(nh3.fc, nh3.kctpm,'-', nh3.fc, nh3.kctum,'-');
  grid on;title('kCARTA NH3');

figure(10);clf;plot(nh3.fc, nh3.kctpm - nh3.kctum,'-');grid on;xlim([650 1800]);
 ylabel('Perturb minus unpert (K)');
 title('kCARTA NH3');

nf11=figure(11);clf;
 h1=subplot(211);hold on; plot(nh3.fc, nh3.kctpm - nh3.kctum,'-');
     plot(frq, srtpm - srtum, '-'); grid on; xlim([750 1150]);
     legend('kCARTA','SARTA');ylabel('pert minus unpert')
     title('SARTA vs kCARTA 10% HN3 col')
 h2=subplot(212);plot(frq, (nh3.kctpm - nh3.kctum) - (srtpm - srtum)','-');
    grid on;xlim([750 1150]);ylabel('kC minus Sarta K');xlabel('wavenumber cm^{-1}')
    linkaxes([h1 h2],'x');set(h1,'xticklabel','');
    pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
    pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])
    
    
%}
