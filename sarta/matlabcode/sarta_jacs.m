% Calculate Jacobians with SARTA

cd /home/chepplew/projects/sarta/cris_hr

addpath /asl/matlib/h4tools                       % rtpread.m
addpath /home/sergio/MATLABCODE                   % replicate_rtp_headprof.m
addpath /asl/matlib/rtptools                      % subset_rtp.m
addpath /asl/matlib/aslutil                       % rad2bt.m

% -------------------------
% 1 - Profile Jacobian
% --------------------------
% Get a regression profile
%fnrtp = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1100_400ppm.op.rtp';
fnrtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_2235.op.rtp';
[hdr har pdr par] = rtpread(fnrtp);
[hnew pnew] = replicate_rtp_headprof(hdr,pdr,49,101);       % 49th std US atm.

% Choose gas
Gases = {'H2O','CO2','O3','N2O','CO','CH4','SO2','HNO3'};
myGas = 'HNO3';

for il = 1 : 100
 %% ioth layer of ith profile will have the perturbation
   switch myGas
   case 'H2O'  
     pnew.gas_1(il,il) = pnew.gas_1(il,il) * 1.1;
   case 'CO2'  
     pnew.gas_2(il,il) = pnew.gas_2(il,il) * 1.1;
   case 'O3'   
     pnew.gas_3(il,il) = pnew.gas_3(il,il) * 1.1;
   case 'N2O'  
     pnew.gas_4(il,il) = pnew.gas_4(il,il) * 1.1;
   case 'CO'   
     pnew.gas_5(il,il) = pnew.gas_5(il,il) * 1.1;
   case 'CH4'  
     pnew.gas_6(il,il) = pnew.gas_6(il,il) * 1.1;
   case 'SO2'  
     pnew.gas_9(il,il) = pnew.gas_9(il,il) * 1.1;
   case 'HNO3' 
     pnew.gas_12(il,il) = pnew.gas_12(il,il) * 1.1;
   end
end
%% 101st profile is unperturbed, so it is the "standard"

% rtpwrite:
dsav  = '/asl/s1/chepplew/projects/sarta/cris_hr/';
fortp = ['jac_' myGas '_r49_pert_1100mb_400pp.rtp'];
rtpwrite([dsav fortp], hnew, har, pnew, par);

% then run sarta, 
sarta_exe = '/home/chepplew/gitLib/sarta/bin/sarta_kc_cris_hrg4_400p_full_th2';
rtp_in    = [dsav fortp];
rtp_out   = [dsav 'jac_' myGas '.rtp'];

if(exist('sarta.out')) delete('sarta.out'); end
unix([ sarta_exe ' fin=' rtp_in ' fout=' rtp_out ' >sarta.out']);

% then re-read
[hdx hax pdx pax] = rtpread(rtp_out);
btx = rad2bt(hdx.vchan, pdx.rcalc);

% then for the channel you want 
switch myGas
  case 'H2O'
  case 'CO2'
  case 'O3'
    ic = find(hdx.vchan > 1040,1) - 1;            % O3: 1040 cm-1
  case 'N2O'
  case 'CO'
  case 'CH4'
    ic = find(hdx.vchan > 1305,1) - 1;            % CH4:1305 cm-1
  case 'SO2'
    ic = find(hdx.vchan > 1351,1) - 1;
  case 'HNO3'
    ic = find(hdx.vchan > 879,1) - 1;
end

% just plot(tc(iC,101)-tc(iC,:)))
figure(1);clf;semilogy( (btx(ic,:)-btx(ic,101))/0.1,pdx.plevs(:,1),'.-' );
 set(gca,'YDir','Reverse');grid on;
 ylim([0.005 1100]);title(['Jacobian ' myGas ' 1350wn']);xlabel('d(BT)/d(SO2)');
 ylabel('Pressure Layer hPa');
% saveas(gcf,['./figs/jac_' myGas '_1350wn_sarta.png'],'png')

% Compare with kcarta results
dpk = '/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/';
KJ = load([dpk 'g6_jac.mat']);
KJ = load([dpk 'g3_jac.mat']);
%KJ = load([dpk 'g12_jac.mat']);
ick = find(KJ.fout>1305,1) - 1;               % CH4
ick = find(KJ.fout>1040,1);               % O3

figure(2);clf;semilogy(KJ.jout(ick,:),pdx.plevs(97:-1:1,1),'.-');grid on;
title('Jac kcarta O3 1040 wn');
set(gca,'YDir','Reverse');ylim([0.005 1100]);
% saveas(gcf,'./figs/jac_O3_1040wn_kcarta.png','png');

% -------------------------------
% 2 - Spectral Jacobian
% -------------------------------
% Get the test perturbed set
fnrtp = '/asl/s1/chepplew/projects/sarta/cris_hr/testperturb_2235.rp.rtp';
[hdr har pdr par] = rtpread(fnrtp);

%
  pind =  {{'unp',  0,  10,20,30,40,50,60},...
           {'WV',   1,  01,11,21,31,41,51},... 
           {'CO2',  2,  02,12,22,32,42,52},...
           {'O3',   3,  03,13,23,33,43,53},...
           {'N2O',  4,  04,14,24,34,44,54},...
           {'CO',   5,  05,15,25,35,45,55},...
           {'CH4',  6,  06,16,26,36,46,56},...
           {'SO2',  9,  07,17,27,37,47,57},... 
           {'HNO3',12,  08,18,28,38,48,58},...
           {'T',   99,  09,19,29,39,49,59}};

uindx = [10 20 30 40 50 60];
% Choose gas:
myGas = 'O3'; mygid = 3;
for ip=1:10 if(ismember(myGas,pind{ip}{1})); ilut=ip; end; end
for i=1:6 gindx(i) = int32(pind{ilut}{i+2}); end;

% then run sarta, 
rtp_in    = [dsav fnrtp];
rtp_out   = [dsav 'sar_testperturb.rtp'];

if(exist('sarta.out')) delete('sarta.out'); end
unix([ sarta_exe ' fin=' rtp_in ' fout=' rtp_out ' >sarta.out']);

% then re-read (already computed)

fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_tp60_kc_400p_full.rtp';
[hds has pds pas] = rtpread(fnsar);

bt_gas = real(rad2bt(hds.vchan,pds.rcalc(:,gindx)));
bt_unp = real(rad2bt(hds.vchan,pds.rcalc(:,uindx)));
jac_gas = (nanmean(bt_gas,2) - nanmean(bt_unp,2))/0.1;

figure(3);clf;plot(hds.vchan,jac_gas,'-');xlim([640 1150]);grid on;
title([myGas ' Jacobian testpurb.rtp']);
xlabel('wavenumber cm^{-1}'); ylabel('Jac d(BT)/d(SO2) K/pc');
% saveas(gcf,['./figs/jac_' myGas '_sarta_spectrum.png'],'png');

