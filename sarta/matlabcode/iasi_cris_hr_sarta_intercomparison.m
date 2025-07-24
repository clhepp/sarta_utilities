% script iasi_cris_hr_sarta_intercomparison()

% Initialize

addpath /asl/matlib/h4tools               % rtpread
addpath /asl/matlib/rtptools              % rtpwrite_12
addpath /asl/matlib/aslutil               % mktemp
addpath /home/motteler/shome/iasi_decon   % iasi2cris
addpath /asl/packages/ccast/source        % inst_params
addpath /asl/packages/airs_decon/source   % hamm_app

% the sarta executables:
sarta.iasi    = '/home/chepplew/gitLib/sarta/bin/iasi_jun19';
sarta.cris_hr = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_prod';

% The test profiles
srcrtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_unitemis_seaemis_7angs_night.rtp';
[head, hatt, prof, patt] = rtpread(srcrtp);
if(isfield(prof,'rcalc'))
  prof=rmfield(prof,'rcalc');
end

subs = struct;
subs.satzen = unique(prof.satzen);
subs.uemis  = sort(unique(prof.emis(9,:)));
jlen = length(prof.satzen)
% switch 686 = 49prf*7ang*2sfc or 2352 = 49prf*8ang*6sfc
switch jlen
  case 686
  for i=1:7
    idx.a{i} =  find(prof.satzen == subs.satzen(i) & prof.emis(9,:) == 1);
    idx.b{i} =  find(prof.satzen == subs.satzen(i) & prof.emis(9,:) < 1);
  end
  subs.names = {'scan angles x 7','2 x surface emissivity: unit and sea'};
end

% Preparation and calcs

% IASI
    x      = load('/home/chepplew/myLib/data/f_iasi.mat');
    freq   = x.f_iasi;
    idchan = x.ichan_iasi;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;        % check the kCARTA prediction surface pressure
    %szz = size(prof.zobs);
    %prof.zobs = 815000.0*ones(1,szz(2));

    tmp = mktemp();
    outfiles = rtpwrite_12(tmp,head,hatt,prof,patt);

    ifn_1 = outfiles{1};     ifn_2 = outfiles{2};
    ofn_3 = [tmp '.sar_1'];  ofn_4 = [tmp '.sar_2'];
    tic
      eval(['! ' sarta.iasi ' fin=' ifn_1 ' fout=' ofn_3 ' > /home/chepplew/logs/sarta/sar_out.log']);
      eval(['! ' sarta.iasi ' fin=' ifn_2 ' fout=' ofn_4 ' > /dev/null']);
    toc
    cfin = [tmp '.sar'];
    [hds,~,ptemp,~] = rtpread_12(cfin);
    calc.irad  = ptemp.rcalc;
    calc.ibsc  = rad2bt(hds.vchan, ptemp.rcalc);
    calc.ifrq  = hds.vchan;
    delete(tmp)
    
% CrIS.FSR
    x       = load('/home/chepplew/myLib/data/fcris_hires_4grd.mat');
    freq    = x.vchan;
    idchan  = x.ichan;
    head.vchan = freq;
    head.ichan = idchan;
    head.nchan = length(idchan);
    %head.pmax  = 1013;
    tmp = mktemp();
    rtpwrite(tmp,head,hatt,prof,patt);

    ifn = tmp;      ofn = [tmp '.sar'];
    tic
      eval(['! ' sarta.cris_hr ' fin=' ifn '  fout=' ofn ' > /home/chepplew/logs/sarta/sar_out.log']);
      %command=[SARTAEXE ' fin=rtpx fout=' sarout ' > /home/chepplew/logs/sarta/sar_out.log'];
      %system(command);
    toc
    [hds,~,ptemp,~] = rtpread(ofn);
    calc.crad = ptemp.rcalc;
    calc.cbsc = rad2bt(hds.vchan, ptemp.rcalc);
    calc.cfrq = hds.vchan;

% Translation IASI to CrIS.FSR
opt1 = struct;
opt1.hapod    = 1;
opt1.ngaurd   = 4;
opt1.user_res = 'hires';

[i2crad, i2cfrq] = iasi2cris(calc.irad, calc.ifrq, opt1);

[i2, i4] = seq_match(i2cfrq, calc.cfrq);


%{
% Plotting

% use first 6 view angles for each surface subset:
iwnt = [];
for j=1:6
  iwnt = [iwnt idx.b{j}];
end
iwnt = sort(iwnt);


figure; hold on;
  plot(calc.ifrq, rad2bt(calc.ifrq, nanmean(calc.irad,2)),'-');
  plot(calc.cfrq, rad2bt(calc.cfrq, nanmean(calc.crad,2)),'-');


figure; hold on;
  plot(i2cfrq(i2), rad2bt(cfreq(i2), nanmean(crad(i2,:),2)),'-');
  plot(calc.cfrq(i4), rad2bt(calc.cfrq(i4), nanmean(calc.crad(i4,:),2)),'-')

bt2m = real(rad2bt(i2cfrq(i2), nanmean(i2crad(i2,iwnt),2)));
bt4m = rad2bt(calc.cfrq(i4),  nanmean(calc.crad(i4,iwnt),2));
    
figure;clf; 
  plot(i2cfrq(i2), bt2m - bt4m,'-')  


%}

%{
% BUILD NOTES:

target->  iasi_jun19:   exec: ../bin/iasi_jun19
include: incFTC_iasi_jun19.f  -> /asl/data/sarta_coef/Data_IASI_jun19/

test: turn off optran
sarta.iasi  = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_notra';

Similarly build CrIS.FSR with no optran:
include: incFTC_cris_hrg4_p2019_dec2018_prod.f
sarta.cris_hr = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18_notra';


%}
