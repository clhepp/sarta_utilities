addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools

fip = 'junk.ip.rtp';
fop = 'junk.op.rtp';
frp = 'junk.rp.rtp';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

newsarta = ['/home/chepplew/gitLib/sarta/bin/sarta_cloudy_airs_may19'];
oldsarta = ['/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003'];
oldsarta = ['/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3'];

klayerser  = ['!' klayers ' fin=' fip ' fout=' fop];
sartaer1   = ['!' oldsarta   ' fin=' fop ' fout=' frp];
sartaer2   = ['!' newsarta   ' fin=' fop ' fout=' frp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /home/sergio/PCRTM_XIANGLEI/RUN2019/SergioOutput/V2/savedata10.31.215.mat
h = savedata.h2x;
p = savedata.p2x;

rtpwrite(fip,h,[],p,[]);
eval(klayerser)

eval(sartaer1)
[h1,ha1,p1,pa1] = rtpread(frp);

eval(sartaer2)
[h2,ha2,p2,pa2] = rtpread(frp);
rmer = ['!/bin/rm junk.ip.rtp junk.op.rtp junk.rp.rtp']; eval(rmer);

tobs = rad2bt(h.vchan,p.robs1);
tcal1 = rad2bt(h.vchan,p1.rcalc);
tcal2 = rad2bt(h.vchan,p2.rcalc);

plot(h.vchan,tcal1-tcal2)
plot(h.vchan,tobs-tcal1,'bx',h.vchan,tobs-tcal2,'ro')

plot(p.cngwat + p.cngwat2,tcal1(5,:),'b.',p.cngwat+p.cngwat2,tcal2(5,:),'r.')
hl = legend('BT1231 old','BT1231 new','location','best'); grid; xlabel('ice+water loading'); ylabel('BT1231')
axis([0 200 250 300])

plot(p.stemp,tcal1(5,:),'b.',p.stemp,tcal2(5,:),'r.',p.stemp,p.stemp,'k')
hl = legend('BT1231 old','BT1231 new','location','best'); grid; xlabel('Stemp'); ylabel('BT1231')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to do clear cals'); pause
p.cngwat = 0*p.cngwat;
p.cfrac  = 0*p.cfrac;
p.ctype  = ones(size(p.ctype)) * -9999;
p.cngwat2 = 0*p.cngwat2;
p.cfrac2  = 0*p.cfrac2;
p.ctype2  = ones(size(p.ctype2)) * -9999;
p.cfrac12  = 0*p.cfrac12;

rtpwrite(fip,h,[],p,[]);
eval(klayerser)

eval(sartaer1)
[h1,ha1,p1,pa1] = rtpread(frp);

eval(sartaer2)
[h2,ha2,p2,pa2] = rtpread(frp);

rmer = ['!/bin/rm junk.ip.rtp junk.op.rtp junk.rp.rtp']; eval(rmer);

tobs = rad2bt(h.vchan,p.robs1);
tcal1 = rad2bt(h.vchan,p1.rcalc);
tcal2 = rad2bt(h.vchan,p2.rcalc);

plot(h.vchan,tcal1-tcal2)
plot(h.vchan,tobs-tcal1,'bx',h.vchan,tobs-tcal2,'ro')

%plot(p.cngwat + p.cngwat2,tcal1(5,:),'b.',p.cngwat+p.cngwat2,tcal2(5,:),'r.')
%hl = legend('BT1231 old','BT1231 new','location','best'); grid; xlabel('ice+water loading'); ylabel('BT1231')
%axis([0 200 250 300])

plot(p.stemp,tcal1(5,:),'b.',p.stemp,tcal2(5,:),'r.',p.stemp,p.stemp,'k')
hl = legend('BT1231 old','BT1231 new','location','best'); grid; xlabel('Stemp'); ylabel('BT1231')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin3 = 'junk3a.rtp';
fop3 = 'junk3b.rtp';
frp3 = 'junk3c.rtp';
klayerser  = ['!' klayers ' fin=' fin3 ' fout=' fop3];


fnrtp = '/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/10/31/cloudy_airs_l1c_ecm_sarta_baum_ice.2018.10.31.222.rtp';

[head, hatt, prof, patt] = rtpread(fnrtp);
head.ptype=0;
rtpwrite(fin3, head,[],prof,[]);
eval(klayerser)

ich = find(head.vchan > 1231,1);

figure(3);clf;
  simplemap(prof.rlat, prof.rlon, rad2bt(head.vchan(ich),prof.rcalc(ich,:)))

sartaer3   = ['!' newsarta   ' fin=' fop3 ' fout=' frp3];
eval(sartaer3)
[h1,ha1,p1,pa1] = rtpread(frp);

