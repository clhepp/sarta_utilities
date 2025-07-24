
cd /home/chepplew/data/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

fnsrc  = '/asl/rtp/rtp_airicrad_v6/clear/2019/era_airicrad_day183_clear.rtp';
fortp  = 'junk_op.rtp';
fortp2 = 'junk_op2.rtp';
fnrtp  = 'junk_ip.rtp';

%hd.ptype=0;

KLAYERSEXE  = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
klayers_run = [KLAYERSEXE ' fin=' fnsrc ' fout=' fnrtp ' > ' ...
               '/home/chepplew/logs/klayers/klout.txt'];
unix(klayers_run);

[hd ha pd pa] = rtpread(fnrtp);

hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
freq     = hdfread(hinfo.SDS(2));
idchan   = hdfread(hinfo.SDS(1));
hd.vchan = single(freq);
hd.ichan = int32(idchan);
hd.nchan = length(idchan);

rtpwrite(fnrtp, hd, ha, pd, pa);


SARTAEXE = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
SARTANEW = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod_tra';

eval(['! ' SARTAEXE ' fin=' fnrtp '  fout=' fortp ' > /dev/null']);
eval(['! ' SARTANEW ' fin=' fnrtp '  fout=' fortp2 ' > /dev/null']);

[~,~,ptemp1,~] = rtpread(fortp);
[~,~,ptemp2,~] = rtpread(fortp2);

bc.old = rad2bt(hd.vchan, ptemp1.rcalc);
bc.new = rad2bt(hd.vchan, ptemp2.rcalc);

figure; plot(hd.vchan, nanmean(bc.old,2),'-')
  hold on; plot(hd.vchan, nanmean(bc.new,2),'-')
clf; plot(hd.vchan, nanmean(bc.old,2) - nanmean(bc.new,2),'-')


