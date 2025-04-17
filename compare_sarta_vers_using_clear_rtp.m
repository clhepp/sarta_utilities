% compare_sarta_vers_using_clear_rtp.m

% For comparing results for a whole day of clear RTP subset using two versions of
% sarta. Used to chack consistency and change of spectroscopy &c.

% 0. Choose which sensor
% 1. Choose which two sarta execs
% 2. Choose which day of clear RTP to use

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /asl/matlib/plotutils                   % aslprint
addpath /asl/packages/airs_decon/source         % hamm_app

% Choose sensors
csens = 'AIRS_L1C';

% KLAYERS 
klayersexe = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';


% Get spectral grid to update head structure and the two sarta execs
switch csens
  case 'AIRS_L1C'
    fn_rtp1 = '/asl/rtp/airs/airs_l1c_v674/clear/2021/ecmwf_airicrad_day270_clear.rtp';
    vchan = hdfread('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf','freq');
    ichan = hdfread('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf','chanid');
    sartaexe1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod_v2';
    sartaexe2 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev';
  case 'CRIS_FSR'
    fn_rtp1 = '/asl/rtp/cris/j01_ccast_hires/clear/2019/cris2_era_csarta_clear_d20190621.rtp';

  case 'CRIS_NSR'

  case 'CHIRP'

  case 'IASI'

end

% Get the original clear subset RTP file:
[head hatt prof patt] = rtpread(fn_rtp1);
f2645 = head.vchan;

% screen for NaNs in emis for landfrac != 0.00 (bug in rtp_add_emis)
iiland = find(prof.landfrac > 0);
for ip = iiland
  iinan = find(isnan(prof.emis(:,ip)));
  while(~isempty(iinan))
    prof.emis(iinan,ip) = prof.emis(iinan+1,ip);
    iinan      = find(isnan(prof.emis(:,ip)));
  end
  iinan = find(isnan(prof.rho(:,ip)));
  while(~isempty(iinan))
    prof.rho(iinan,ip) = prof.rho(iinan+1,ip);
    iinan      = find(isnan(prof.rho(:,ip)));
  end
end

% Update head struct and create temporary files
head.vchan = vchan;
head.ichan = ichan;
head.nchan = length(ichan);

% and save ready for processing
tmp = mktemp();

fn_1 = [tmp '.ip.rtp'];
fn_2 = [tmp '.op.rtp'];
fn_3 = [tmp '.sar.rtp'];

rtpwrite(fn_1,head,hatt,prof,patt);

klayers_cmd = [klayersexe ' fin=' fn_1 ' fout=' fn_2 ' > ' ...
               '/home/chepplew/logs/klayers/klout.txt'];
system(klayers_cmd)

SARTAEXE = sartaexe1;
sarta_cmd = [SARTAEXE ' fin=' fn_2 ' fout=' fn_3 ' > /home/chepplew/logs/sarta/sarta_out.txt'];
system(sarta_cmd)
[hds1, has1 ,pds1 ,pas1] = rtpread(fn_3);

% Repeat calc using second sarta exec
SARTAEXE = sartaexe2;
sarta_cmd = [SARTAEXE ' fin=' fn_2 ' fout=' fn_3 ' > /home/chepplew/logs/sarta/sarta_out.txt'];
system(sarta_cmd)
[hds2, has2 ,pds2 ,pas2] = rtpread(fn_3);

% Compute BT spectra
[sfreq,iss] = sort(hds1.vchan);
[ix iy]     = seq_match(f2645, sfreq);

bto   = rad2bt(f2645, prof.robs1);
btsc1 = rad2bt(hds1.vchan, pds1.rcalc);
btsc2 = rad2bt(hds2.vchan, pds2.rcalc);

iiwnt = find(pds1.solzen > 90);
plot(sfreq, rad2bt(sfreq, mean(pds1.rcalc(iss,iiwnt),2))- ...
     rad2bt(sfreq, mean(pds2.rcalc(iss,iiwnt)>> ),'-')

plot(f2645, rad2bt(f2645,mean(prof.robs1(:,iiwnt),2)),'-')
plot(sfreq(iy), rad2bt(sfreq(iy),mean(pds1.rcalc(iy,iiwnt),2)),'-')


%{
% Plotting
% - maps


% - spectra

bands = [650, 1090; 1200, 1745; 2150, 2560];
figure(1);clf;plot(fcris1, bom, '-', fcris1, bc1m,'-',fcris2,bc2m,'-' );
  axis([650 1745 215 280]);grid on;
  title('CrIS hiRes clear sub 2016d020 Obs vs calc');xlabel('wn cm-1');ylabel('BT (K)');
  legend('Obs','old sarta','new sarta')
figure(2);clf;plot(fcris1, bo2m,'-', fcris2, bc2m,'-' );axis([650 1745 215 280]);grid on;
  title('Obs (b)  calc (r) LW new SARTA');xlabel('wn cm-1');ylabel('BT (K)');
 % aslprint('./figs/cris_d2016020_clr_sarta_btMean_LW_MW.png');

for ib=1:3
figure(ib);clf;
  h1=subplot(2,1,1);plot(freq,bc2m,'-',freq,bo2m,'-');grid on;xlim([bands(ib,:)]);
  h2=subplot(2,1,2);plot(freq, bo2m - bc2m,'-' );axis([bands(ib,:) -2 2]);grid on;
  title('CrIS hiRes clr sub 2016d020 bias: obs - calc');xlabel('wn cm-1');ylabel('BT (K)');
  hold on;plot(fcris1(1:end-12), bom(1:end-12) - bc2m(13:end),'-' );
  legend('old sarta','new sarta');
% LW: fcris1(1:end-4), bc2m(5:end). MW: fcris1(1:end-12), bc2m(13:end);
 % aslprint('./figs/cris_2016d020_clr_sarta_bias_MW.png');
end
%}
