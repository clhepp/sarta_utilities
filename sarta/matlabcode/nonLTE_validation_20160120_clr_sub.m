 addpath /asl/matlib/plotutils
 addpath /asl/packages/airs_decon/source

 [hds has pds pas] = rtpread('/asl/rtp/rtp_cris_ccast_hires/clear/2016/cris_hr_ecmwf_d20160120_clear_csarta.rtp');
>> [nx ny] = size(pds.plevs)

inx = find(pds.landfrac == 0 & abs(pds.rlat) < 40.0);
ind = find(pds.landfrac == 0 & abs(pds.rlat) < 40.0 & pds.solzen < 90);
inn = find(pds.landfrac == 0 & abs(pds.rlat) < 40.0 & pds.solzen >= 90);

figure(1);clf;simplemap(pds.rlat(ind),pds.rlon(ind),pds.robs1(472,ind));

btcd = real(rad2bt(hds.vchan,pds.rcalc(:,ind)));
btcn = real(rad2bt(hds.vchan,pds.rcalc(:,inn)));
btcdm = nanmean(btcd,2);
btcnm = nanmean(btcn,2);


  junk = single(hamm_app(double(pds.robs1(:,ind))));
btod = real(rad2bt(hds.vchan,junk));
  junk = single(hamm_app(double(pds.robs1(:,inn))));
bton = real(rad2bt(hds.vchan,junk));
btodm = nanmean(btod,2);
btonm = nanmean(bton,2);

dbias = btodm - btcdm;
nbias = btonm - btcnm;
dstdv = nanstd(btod - btcd,1,2);
nstdv = nanstd(bton - btcn,1,2);

figure(2);clf;plot(hds.vchan,btcdm,'-',hds.vchan,btcnm,'-');grid on; 
   legend('day','night');xlim([2150 2550]);
   title('csarta calc, tropical ocean scenes, day/night');

figure(3);clf;plot(hds.vchan,btodm,'-',hds.vchan,btonm,'-');grid on; 
   legend('day','night');xlim([2150 2550]);
   title('csarta obs, tropical ocean scenes, day/night');

figure(4);clf;plot(hds.vchan,dbias,'-',hds.vchan,nbias,'-');grid on;hold on;
  plot(hds.vchan,dstdv,'-',hds.vchan,nstdv,'-');
  legend('day bias','night bias','day stdv','night stdv');xlim([2155 2550]);
  title('mean bias and stddev obs-calc trop ocean');xlabel('wn cm^{-1}');ylabel('BT K');  
  % saveas(gcf,'./figs/d20160120_tropOcean_clr_day_night_bias_stdv_SW.png','png');
