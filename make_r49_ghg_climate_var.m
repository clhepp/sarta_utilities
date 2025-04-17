% make_r49_ghg_climate_var.m
%
% duplicates original 49 regression profiles for fast coefficient regression
%   with CO2, CH4, N2O varying abundance to simulate 20+ (40) years of
%   climate growth.
%   CH4 1775:1925 ppbb
%
%
%
%
%
%

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt, int2bits, mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/chepplew/projects/sarta/matlabcode
addpath /home/cheppelw/myLib/matlib/rtptools     % replicate_head_prof
aadpath /home/motteler/shome/airs_decon/source   % seq_match

% Avogadros no.
AVO = AVO=6.02214076E23;

% No. density of air at STP: (molec.m-3)
Na = 2.69E25;

% RMM CO2, CH4, N2O
ch4_rmm = 16.0425; % g/mol
co2_rmm = 44.0095; % g/mol
n2o_rmm = 44.0128; % g/mol

% KLAYERS
klayers_exe = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';

% the five CO2, CH4 and N2O test-point values:
co2ppm = [360:15:420];
ch4ppb = [1775:37.5:1925];
n2oppb = [314:5.5:336];

% -------------------------------------------------------------------------------
% Original R49 RTP:
% Defaults to: satzen=0, solzen=150. emis=sea emissivity. spres = 1013.25
srcfn1 = '/asl/s2/hannon/SpinachHome/Fit_deltaR_nonLTE/pin_feb2002_reg_op.rtp';

srcdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
srcfn2 = [srcdr  'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];
srcfn3 = [srcdr 'REGR49_400ppm_H2012_June2016/regr49_1013_400ppm.ip.rtp'];

% Paths too long for input param to klayers.
system(['cp ' srcfn1 ' .'])

[head hatt prof patt] = rtpread(srcfn1);
[nr np] = size(prof.emis);

% srcfn1: spres = 1100mb nlevs=101 used for OD regression w/44 gases
% Need to add head.{v,i,n}chan
[hdX,~,pdX,~] = rtpread(srcfn2);
head.ichan = hdX.ichan;
head.vchan = hdX.vchan;
head.nchan = hdX.nchan;

% Fill gases CH4 N2O
rdate = '2012/01/01';
tdate = dnum2tai(datenum(rdate,'yyyy/mm/dd'))
ptime = tdate*ones(size(prof.pobs));
prof.ptime = ptime;
prof.rtime = ptime;
prof.rlat = ones(size(prof.pobs));
prof.rlon = ones(size(prof.pobs));
[head,hatt,prof] = fill_co2(head,hatt,prof);
[head,hatt,prof] = fill_n2o(head,hatt,prof);
[head,hatt,prof] = fill_ch4(head,hatt,prof);
%
% this input rtp already has gas_6 (CH4) so need to correct head.glist
head.glist = head.glist([1:3 5:8]);
head.gunit = head.gunit([1:3 5:8]);
head.ngas  = 7;


% For N2O need to fill missing profile data:
% use klayers to fill with AFGL N2O
tmpPath = mktemp();
fn_ip   = srcfn;
fn_op  = [tmpPath '_kl.op.rtp'];
command = [klayers_exe ' fin=' fn_ip ' fout=' fn_op ' > /home/chepplew/logs/klayers/klout.txt'];
system(command)
%
[hdk,~,pdk,~] = rtpread(fn_op);



% Create duplicates centered around [360:15:420] co2ppm
% use klayers to populate prof.gas_2 w/ co2ppm values prescribed.

% Calc fractional changes. Use SI meters..
% get layer thickness, use sfc thickness
thick = mean(prof.palts(1:end-1,:)-prof.palts(2:end,:),2);
% %ch4frac = ch4ppb./mean(prof.gas_6(1,:),2);
% %n2ofrac = 1E-3*n2oppb./mean(prof.gas_4(1,:),2);
co2frac = 1E-6*thick(end)*co2ppm*Na./(1E4*mean(prof.gas_2(98,:),2));
ch4frac = 1E-9*thick(end)*ch4ppb*Na./(1E4*mean(prof.gas_6(98,:),2));
n2ofrac = 1E-9*thick(end)*n2oppb*Na./(1E4*mean(prof.gas_4(98,:),2));

h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,5);
   % %py.gas_2(:,1)  = 1.0   * prof.gas_2(:,ii);          % CO2
   py.co2ppm(1:5)   = co2ppm; 
   py.gas_2(:,1:5)  = co2frac.*prof.gas_2(:,ii);          % CO2
   py.gas_6(:,1:5)  = ch4frac.*prof.gas_6(:,ii);          % CH4
   py.gas_4(:,1:5)  = n2ofrac.*prof.gas_4(:,ii);          % N2O
   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2, hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n');

% trim out data not needed (and kCARTA complains)
[hdo,pdo] = subset_rtp(h2,p2,[1 2 3 4 5 6 9 10 11 12 20],[],[]);

% RTA now needs prof.udef(20,:) set
udef_20 = zeros(size(pdo.nlevs));
pdo.udef(20,:) = udef_20;

% RTA needs emissivity data
pdo.efreq   = repmat(pdX.efreq,1,5);
pdo.emis    = repmat(pdX.emis,1,5);
pdo.rho     = repmat(pdX.rho,1,5);
pdo.nemis   = repmat(pdX.nemis,1,5);
pdo.satzen  = repmat(pdX.satzen,1,5);
pdo.solzen  = repmat(pdX.solzen,1,5);
pdo.zobs    = repmat(pdX.zobs,1,5);
pdo.scanang = repmat(pdX.scanang,1,5);
pdo.stemp   = repmat(pdX.stemp,1,5);

rtpwrite(fn_op2, hdo, hatt, pdo, patt); 
[hds,~,pds,~] = rtpread(fn_op2);

[sf2378,isy] = sort(hds.vchan);
sar.btc = rad2bt(sf2378, pds.rcalc(isy,:));

% ================================================
/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev
fn_out = [tmpPath '_sar_out.rtp'];
sar_cmd = [SARTAEXE ' fin=' fn_op2 ' fout=' fn_out ' > /home/chepplew/logs/sarta/sar_out.txt'];

system(sar_cmd)
[hds,has,pds,pas] = rtpread(fn_out);

% ===========================================================

% Update CO2 gas amount profiles by running through klayers with upd. co2ppm
fn_ip2 = [tmpPath '_kl2.ip.rtp'];
fn_op2 = [tmpPath '_kl2.op.rtp'];
rtpwrite(fn_ip2, h2, hatt, p2, patt); 

command = [klayers_exe ' fin=' fn_ip2 ' fout=' fn_op2 ' > /home/chepplew/logs/klayers/klout.txt'];
system(command)
%
[hdk,~,pdk,~] = rtpread(fn_op2);

% Record Index of each series
ns = length(p2.pobs);
idx = struct
idx(1).sc = [1:5:ns];
idx(2).sc = [2:5:ns];
idx(3).sc = [3:5:ns];
idx(4).sc = [4:5:ns];
idx(5).sc = [5:5:ns];






%{
fnin  = '/home/chepplew/projects/sarta/regr49_1013_400ppm.ip.rtp';
fnout = './r49.op.rtp';
command = [klayers_exe ' fin=' srcfn ' fout=' fnout ' > /home/chepplew/logs/klayers/klout.log'];


%}

kc.home = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/';
kc.home = '/home/chepplew/data/sarta/prod_2023/generic/kcarta/245_ghg_tests/';
kc.list = dir([kc.home 'individual_prof_convolved_kcarta_airs*']);
for ifn=1:245
  xx=load([kc.list(ifn).folder '/' kc.list(ifn).name]);
  junk = strsplit(kc.list(ifn).name,{'_','.'});
  pnum(ifn) = str2num(junk{6});
  kc.rad(:,ifn) = xx.rKc;
end
[sfreq,iss] = sort(xx.fKc);
kc.btc = rad2bt(sfreq, kc.rad(iss,:));
% re-order kcarta.dat files as returned by dir()
[spnum ips] = sort(pnum);
kc.btc = kc.btc(:,ips);

% for ease of comparison match kc.bt to sar.bt frequencies
[ix iy] = seq_match(sfreq, sf2378);
kc.btci = kc.btc(ix,:);

% plot(sf2378(iy), mean(sar.btc(iy,:),2),'-', sfreq(ix), mean(kc.btc(ix,:),2),'-')
% plot(sf2378, mean(sar.btc,2),'-', sfreq(ix), mean(kc.btci,2),'-')

% variation with climate - include quadratic fit
X1 = co2ppm - co2ppm(3);
X2 = ch4ppb - ch4ppb(3);
X3 = n2oppb - n2oppb(3);

X = X2;

ichs = find(sf2378>710,5);
sar.btc_mn = []; kc.btc_mn = [];
for i = 1:5
  sar.btc_mn(i) = rms(sar.btc(ichs,idx(i).sc),[1 2]);
  kc.btc_mn(i)  = rms(kc.btci(ichs,idx(i).sc),[1 2]);
end

sar.btc_mn = []; kc.btc_mn = [];
for ich = 1:length(sf2378)
  for i = 1:5
    sar.btc_mn(i,ich) = mean(sar.btc(ich,idx(i).sc), 2);
    kc.btc_mn(i,ich)  = mean(kc.btci(ich,idx(i).sc), 2);
  end
end

Y = []; sar.coefs = []; kc.coefs = [];
for ich = 1:length(sf2378)
  Y      = sar.btc_mn(:,ich)-mean(sar.btc_mn(:,ich),1);
  fmodel = fit(X(:), Y(:),'poly2');
  sar.coefs(ich,:) = [fmodel.p1 fmodel.p2 fmodel.p3];
  Y      = kc.btc_mn(:,ich)-mean(kc.btc_mn(:,ich),1);
  fmodel = fit(X(:), Y(:),'poly2');
  kc.coefs(ich,:) = [fmodel.p1 fmodel.p2 fmodel.p3];
end

% plot(sf2378, sar.coefs(:,2),'b-', sf2378,sar.coefs(:,1),'r-'); hold on;
% plot(sf2378, kc.coefs(:,2),'c-',  sf2378,kc.coefs(:,1),'m-')
% legend('sar.lin','sar.quad*10','kc.lin','kc.quad*10','location','southEast')
% title('quad fit coefs <gas> var SARTA. kCARTA')
% xlabel('wavenumber cm-1');ylabel('poly fit coefficient (/K)')

ichs = find(sf2378>1300,5);

figure; plot(mean(sar.btc(ichs,id


%}
