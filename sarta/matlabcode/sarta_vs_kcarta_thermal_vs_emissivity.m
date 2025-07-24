% sarta_vs_kcarta_thermal_vs_emissivity

% Run sarta twice with new exe and airs.apr08 version
% structure: calc =
%    bsc_ch: [2834x2352 single]
%    bsc_sc: [2834x2352 single]
%       bkc: [2834x2352 single]
%       frq: [2834x1 single]
% Select params for LW window channel
ic=find(head.vchan(ib) > 960.0,1)
% a water line in the window region:
ic=find(head.vchan(ib) > 849.4,1);

iwant = [indx.emis1(1) indx.emis2(1) indx.emis3(1) indx.emis4(1) ...
         indx.emis5(1) indx.emis6(1)]
ewant = [prof.emis(6,indx.emis1(1)) prof.emis(6,indx.emis2(1)) ...
   prof.emis(1,indx.emis3(1)) prof.emis(1,indx.emis4(1)) ...
   prof.emis(6,indx.emis5(1)) prof.emis(1,indx.emis6(1))]




figure(4);clf;plot(ewant, krc(ib(ic),iwant),'.')
hold on;
  plot(ewant, ptemp.rcalc(ib(ic),iwant),'o')
  plot(ewant, prof1.rcalc(ib(ic),iwant),'d')
figure(4);clf;plot(ewant, calc.bkc(ib(ic),iwant),'.')
hold on;
  plot(ewant, calc.bsc_sc(ib(ic),iwant),'o')
  plot(ewant, calc.bsc_ch(ib(ic),iwant),'d')
  xlim([0.7 1.05]); legend('kcarta','airs apr09','new sarta')
  xlabel('sfc emissivity');ylabel('toa rad');
  title(['TOA rad predict ch' sprintf('%3d %8.3f',ic, head.vchan(ib(ic))) 'wvn'])


  
