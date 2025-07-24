% script compare_old_sarta_notune.m
%
% To see effect of tuning using sample clear rtp.
% AIRS08 Sarta original and rebuilt with no tuning.
%
%
cd /home/chepplew/data/sarta

addpath /asl/matlib/h4tools         % rtpread
addpath /asl/matlib/aslutil         % rad2bt

fnsrc = '/asl/rtp/rtp_airicrad_v672/clear/2019/era_airicrad_day183_clear.rtp';
[hd ha pd pa] = rtpread(fnsrc);

sarta.clear.orig = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
sarta.clear.notune = '/home/chepplew/projects/sarta/v108_copy/Bin/sarta_apr08_m140_wcon_nte_rbld_notune';

[head hatt prof patt] = rtpread(fnrtp);

% run klayers
fortp       = 'klayers_op.rtp';
KLAYERSEXE  = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
klayers_run = [KLAYERSEXE ' fin=' fnsrc ' fout=' fortp ' > ' ...
               '/home/chepplew/logs/klayers/klout.txt'];
unix(klayers_run);

% inspect
[hd ha pd pa] = rtpread(fortp);


% Run SARTA
  SARTAEXE = sarta.clear.orig;
  ifn = fortp;      ofn = ['sar_op.rtp'];
  eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
  [~,~,ptemp1,~] = rtpread(ofn);
  calc.orig.rad = ptemp1.rcalc;
  calc.orig.bsc = rad2bt(head.vchan, ptemp1.rcalc);
%

  SARTAEXE = sarta.clear.notune;
  ifn = fortp;      ofn = ['sar_op_notune.rtp'];
  eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
  [~,~,ptemp2,~] = rtpread(ofn);
  calc.notune.rad = ptemp2.rcalc;
  calc.notune.bsc = rad2bt(head.vchan, ptemp2.rcalc);
%

% clear unwanted files
delete(fortp)
delete('sar_op.rtp')
delete('sar_op_notune.rtp')

%{
% plots
fa = head.vchan;
figure; plot(fa, rad2bt(fa, nanmean(calc.orig.rad,2)),'-');
  hold on; plot(fa, rad2bt(fa, nanmean(calc.notune.rad,2)),'-')
clf; 
  plot(fa, rad2bt(fa, nanmean(calc.orig.rad,2)) - ...
           rad2bt(fa, nanmean(calc.notune.rad,2)),'-')


%}
