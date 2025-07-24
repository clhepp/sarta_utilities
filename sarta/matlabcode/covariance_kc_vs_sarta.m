function [res] = covariance_kc_vs_sarta(csens);

% covariance_kc_vs_sarta.m
% SYNOPSIS: res = covariance_kc_vs_sarta(csens)
%
% sarta vs kcarta channel covariance for SAF704 profiles
%
% INPUT: csens. Sensor. One of {'AIRS_L1C','CRIS_HR','IASI'};
%
% OUTPUT: res: bias correlation
%
%

addpath /asl/matlib/h4tools                % rtpread
addpath /asl/matlib/rtptools               % rtpwrite_12
addpath /asl/matlib/aslutil                % rad2bt, int2bits, mktemp
addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil                % int2bits
addpath /home/chepplew/projects/sarta/matlabcode

% Colormap for the correlation plts:
load('/home/strow/Matlab/Cmaps/llsmap5.mat');    % colormap(llsmap5);

% Check which sensor
allsens = {'AIRS_L1C','IASI','CRIS_HR'};
csens   = upper(csens);
if(~ismember(csens,allsens)) error('Invalid sensor'); return; end


% All angles (LW use 1:8) and SW)
angles = [0    8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
        65.0428   69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];

% Get the kCARTA TOA Rads
kc.home = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/SAF704_400ppm_H2016_Dec2018_AIRS2834/';
%kc.path = [kc.home 'RAD1013_unitemis/'];
kc.path = [kc.home 'RAD1100_unitemis/'];

kc.dir = dir([kc.path '/convolved_kcarta_RAD*_radiances.mat']);

clear prfnums fnparts;
% reorder these in profile number order
for i=1:length(kc.dir)
  fnparts{i} = strsplit(kc.dir(i).name, {'_' '.'});
  prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
end
[B IB] = sort(prfnums);

disp('Loading kCARTA data')
kc.rad = [];
for fn = IB
  x = load([kc.dir(fn).folder '/' kc.dir(fn).name]);

  switch csens
    case 'AIRS_L1C'
      kc.rad = [kc.rad; x.rairs_all'];

    case 'IASI'
      kc.rad = [kc.rad; x.riasi_all'];
  
    case 'CRIS_HR';
      kc.rad = [kc.rad; x.rcris_all'];
  end
  fprintf(1,'.')
end
fprintf(1,'\n');

% Get spectral grid and SARTA executable:
switch csens
  case 'AIRS_L1C'
    kc.frq = x.fairs;
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_may19';
  case 'IASI'
    kc.frq = x.fiasi;
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19';
  case 'CRIS_HR'
    kc.frq = x.fcris;
    SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16';
end
clear x;
kc.bt = rad2bt(kc.frq, kc.rad');


% Get the SARTA TOA Rads
disp(' Computing SARTA Rads, takes a few minutes')

fnrtp = [kc.home 'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
[head, hatt, prof, patt] = rtpread(fnrtp);

% duplicate profiles at 8 LW angles:
h3 = struct;
p3 = struct;
nprof = size(prof.satzen,2);     % expecting 704
for ii = 1 : nprof
   [hy,py] = replicate_rtp_headprof(head,prof,ii,8);
   for jj = 2:8
     py.satzen(jj) = angles(jj);
   end
  if ii == 1
    h3 = hy;
    p3 = py;
  else
    [h3, p3] = cat_rtp(h3, p3, hy, py);
  end
  %fprintf(1,'.')
  if(~mod(ii,70)); fprintf(1,'.'); end
end
% get rid of rcalc
if(isfield(p3,'rcalc')) p3=rmfield(p3,'rcalc'); end

switch csens
 case 'AIRS_L1C'
    hinfo    = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
    freq     = hdfread(hinfo.SDS(2));
    idchan   = hdfread(hinfo.SDS(1));
    h3.vchan = single(freq);
    h3.ichan = int32(idchan);
    h3.nchan = length(idchan);
      %szz = size(p3.zobs);
      %p3.zobs = 705000.0*ones(1,szz(2));
    tmp = mktemp();
    rtpwrite(tmp,h3,hatt,p3,patt);

    ifn = tmp;      ofn = [tmp '.sar'];
    eval(['! ' SARTAEXE ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);

   [~,~,ptemp,~] = rtpread(ofn);
    sar.frq = h3.vchan;
    sar.rad = ptemp.rcalc;
    sar.bt  = rad2bt(h3.vchan, ptemp.rcalc);

  case 'IASI'
    y        = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq     = y.f_iasi;
    idchan   = y.ichan_iasi;
    h3.vchan = freq;
    h3.ichan = idchan;
    h3.nchan = length(idchan);
      szz = size(p3.zobs);
      p3.zobs = 815000.0*ones(1,szz(2));
    tmp = mktemp();
    outfiles = rtpwrite_12(tmp,h3,hatt,p3,patt);

    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];

    eval(['! ' SARTAEXE ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sarta/sar_out.log']);
    eval(['! ' SARTAEXE ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);

    cfin = [tmp '.sar'];
    [~,~,ptemp,~] = rtpread_12(cfin);
    sar.frq = h3.vchan;
    sar.rad = ptemp.rcalc;
    sar.bt  = rad2bt(h3.vchan, ptemp.rcalc);

  case 'CRIS_HR'
  
end

%
if(strcmp(csens,'AIRS_L1C'))
  [~,ib] = sort(kc.frq);
else
  ib = ':';
end

%{
figure(1);clf;plot(sar.frq(ib),mean(sar.bt(ib,1:8:end),2),'-')
   hold on;plot(kc.frq(ib),mean(kc.bt(ib,1:8:end),2),'-')
%}

disp(' Computing covariance')

kc.cov   = cov(kc.bt');
sar.cov  = cov(sar.bt');
kc.corr  = corrcoef(kc.bt');
sar.corr = corrcoef(sar.bt');
%
% Subset view angles with iang set from 1,2,3,4,5,6 or 7, or 1:7. don't use kCARTA 8th.
iix  = size(prof.satzen,2);
disp([ 'size(prof.satzen,2) ' num2str(iix) ])
iang = [1:7];
oox = [];
for i = 1:length(iang)
  oox = [oox iix(i:8:end)];
end
oox = sort(oox);

bias.cov  = cov( (sar.bt(ib,oox) - kc.bt(ib,oox))' );
bias.corr = corrcoef( (sar.bt(ib,oox) - kc.bt(ib,oox))' );
res = bias.corr;

%{
figure(2);clf;h1=pcolor(kc.frq(ib),kc.frq(ib),bias.corr);shading flat;
   colorbar;caxis([-1.0 1.0]);xlabel('wn cm^{-1}');ylabel('wn cm^{-1}');colormap(llsmap5);
   title('corrcoef SAR minus kc SAF704 at 7 view angles')
%}
