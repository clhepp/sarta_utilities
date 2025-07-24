% compare_kc_sar_r49_prod_2018

% ------------------------------------
% Compare SARTA TOA radiance w/ kcarta
% ------------------------------------
cd /home/chepplew/projects/sarta/prod_2018/

addpath /asl/matlib/h4tools
addpath /asl/packages/airs_decon/source          % hamm_app.m
addpath /asl/matlib/aslutil                      % rad2bt.m
addpath /asl/matlib/plotutils                    % aslprint
addpath /home/chepplew/projects/sarta/matlabcode
warning 'off'

% Choose which regression set to compare.
all_cregr = {'r49','saf704'};
cregr     = 'saf704';

% Enter which sensor to compare w/kcarta {'AIRS' or 'CRIS'}
csens = 'CRIS';

% Default scan angles (LW) for prod_2018 use first 7 (SARTA limited to < 63.deg).
scang = [0   28.6101   38.5184   45.2214   46.3539   48.3250];
scang = [0    8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
        65.0428   69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];
scang = scang(1:7);
for i=1:7 pscan{i} = sprintf('%2.0f',scang(i)); end
iang  =  [1:numel(scang)];

%  1:              BASELINE TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------
% Get the KLAYERS RTP file used to generate the SARTA and kcarta data
% -------------------------------------------------------------------
kpath='/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018/';

%if(strcmp(cregr,'r49'))
  fnlst = dir([kpath 'RAD1013_unitemis/convolved_kcarta_RAD1013_*_radiances.mat']);
% reorder these in profile number order
  for i=1:49 fnparts{i} = strsplit(fnlst(i).name, '_'); 
    prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
  end
  [B IB] = sort(prfnums);  
  krc = []; 
  for i=IB 
    x   = load([kpath 'RAD1013_unitemis/' fnlst(i).name]);
    krc = [krc x.rcris_all];                 % [2235 x 392] 49 profs x 8 angles
  end

% load the SARTA TOA radiances from the r49 profiles (6-angles ea.)
spath   = '/home/chepplew/data/sarta/prod_2018/sarta_data/';
sar_rtp = 'sar_crisg4_basic_mar18_r49_8angs_1013_400_unitemiss.rtp';
sar_rtp = 'sar_crisg4_mar18_basic_optr_co2_r49_7angs_1013_400_unitemiss.rtp';
sar_rtp = 'sar_crisg4_mar18_basic_optr_co2_nte_r49_7angs_1013_400_unitemiss.rtp';
sar_rtp = 'sar_crisg4_mar18_optr_co2_nte_r49_7a_23sz_1013_400_1e.rtp';
sar_rtp = 'sar_crisg4_mar18_basic_so2_r49_7angs_1013_400_1emis.rtp';
sar_rtp = 'sar_crisg4_mar18_basic_so2_7a_1013mb_0p8e.rtp';
[hds has pds pas]     = rtpread([spath sar_rtp]);
[hds2 has2 pds2 pas2] = rtpread([spath sar_rtp]);  % compare two sarta runs.

% load the original 49 profiles input RTP to kcarta and sarta.
r49_rtp = 'regr49_7angs_1013_400ppm_unitemiss.op.rtp';
[hdr har pdr par] = rtpread([spath r49_rtp]);          % 7-angles (restricted < 63-deg)

%      omit ang 8 from kbc
krc_sub = [];
for k=1:49 
  krc_sub = [krc_sub krc(:,[(k-1)*8+1:(k-1)*8+7]) ];
end

% convert to BT and prep for analysis
kbc = rad2bt(x.fcris, krc_sub);
sbc = rad2bt(hds.vchan, pds.rcalc);
sbc2 = rad2bt(hds.vchan, pds2.rcalc);

% calc gobal means for 49 profs & for each of the common 7 scan angles, 
kbcm   = nanmean(kbc,2);         kbcsd = nanstd(kbc,0,2);
sbcm   = nanmean(sbc,2);         sbcsd = nanstd(sbc,0,2);
sbcm2  = nanmean(sbc2,2);       sbcsd2 = nanstd(sbc2,0,2);
radsd  = nanstd(krc_sub - pds.rcalc,0,2);
 ksbm    = 0.5*( kbcm + sbcm );
 mdr     = 1E-3*( 1./drdbt(x.fcris,ksbm) );
btstd    = mdr.*radsd;

% reshape arrays to [nchans x 7 nangs x 49 nprofs] (SARTA is limited to < 63.degs)
kjunk = reshape(kbc,2235,8,[]);
sjunk = reshape(sbc,2235,7,[]); 
kbcm_a = nanmean(kjunk(:,1:7,:),3);
sbcm_a = nanmean(sjunk,3);


%{
% Plot checks
figure(1);clf;plot(x.fcris,kbc(:,1),'-', hds.vchan,sbc(:,1),'-');grid on;
   legend('kcarta','sarta');
nf1=figure(1);clf;set(nf1,'Resize','Off');%set(nf1,'Position',nf1.Position+[0 0 0 420]);
  h1=subplot(411);plot(x.fcris,kbcm,'-', hds.vchan,sbcm,'-');grid on;xlim([640 2600]);
   legend('kcarta','sarta');title('R49 kcarta:sarta.basic.optr.co2 mean spectrum');
  h2=subplot(412);plot(x.fcris,kbcsd,'-', hds.vchan,sbcsd,'-');grid on;xlim([640 2600]);
    legend('kcarta','sarta');title('std.dev spectrum')
  h3=subplot(413);plot(x.fcris,kbcm - sbcm,'-');grid on;legend('kcarta minus sarta');
    title('mean bias');xlim([640 2600]);
  h4=subplot(414);plot(x.fcris, btstd, '-');grid on;title('std.dev of bias (K)');
    xlim([640 2600]);xlabel('wavenumber cm^{-1}');
    
figure(3);clf;hold on; plot(pds.stemp,'bo-'); plot(pdk.stemp,'r+-'); grid on;

figure(3);clf;hold on;plot(sbc(403,1:7:end),'.-');plot(kbc(403,1:8:end),'o-');grid on;
  xlim([1 50]);legend('sarta','kcarta');

figure(4);clf;h1=subplot(211);hold on;grid on;
 for k=1:7 plot(x.fcris, kbcm_a(:,k),'-', hds.vchan,sbcm_a(:,k),'-');end
   legend(pscan,'Location','south','orientation','horizontal');
   xlim([640 1100]);
 h2=subplot(212);hold on;grid on;
 for k=1:7 plot(x.fcris, kbcm_a(:,k) - sbcm_a(:,k),'-');end
   legend(pscan,'Location','south','orientation','horizontal');
   xlim([640 1100]);
   

%}

%{
% Supplement: Compare kcarta HITRAN2012 and HITRAN2016 TOA rads for regr49
kc_path = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';

kc(2).path = [kc_path 'REGR49_400ppm_H2012_June2016/RAD1013/'];
kc(2).list = dir([kc(2).path 'convolved_kcarta_RAD1013_*_radiances.mat']);

kc(3).path = [kc_path 'REGR49_400ppm_H2016_Mar2018/RAD1013_unitemis/'];
kc(3).list = dir([kc(3).path 'convolved_kcarta_RAD1013_*_radiances.mat']);

for m=2:3
  for i=1:49 fnparts{i} = strsplit(kc(m).list(i).name, '_'); 
    prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
  end
  [B IB] = sort(prfnums);  
  kc(m).rc = []; 
  for i=IB 
    x   = load([kc(m).path '/' kc(m).list(i).name]);
    kc(m).rc = [kc(m).rc x.rcris_all];                 % [2235 x 392] 49 profs x 8 angles
  end
  if(m == 2)        % 6 angles
    kc(m).bc = rad2bt(x.fcris, kc(m).rc);
  else if(m == 3)   % 8 angles  reduce to first 6 angles.
    idx = [];                                         
    for k=1:49 idx = [idx [(k-1)*8+1:(k-1)*8+6]]; end  
    kc(m).bc = rad2bt(x.fcris, kc(m).rc(:,idx));
  end
  
  kc(m).bcm = nanmean(kc(m).bc,2); 
end  
%
figure(2);clf;plot(x.fcris, kc(2).bcm,'-', x.fcris,kc(3).bcm,'-');
figure(2);clf;plot(x.fcris, kc(2).bcm - kc(3).bcm,'-');grid on; ylabel('d(BT) K');
 title('RAD1013 H12.r49 minus H16.r49 mean BT (K)');xlabel('wvno. cm^{-1}');

%}
