% make_r49_test_profiles.m
%
%
%

cd /home/chepplew/projects/sarta/prod_2019/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                      % rad2bt, int2bits, mktemp
addpath /asl/matlib/plotutils                    % aslprint
addpath /asl/matlib/rtptools                     % rtpwrite_12
addpath /home/chepplew/projects/sarta/matlabcode


% -------------------------------------------------------------------------------
% Original LEVELS R49 rtp:
srcdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
srcfn = [srcdr 'pin_feb2002_sea_airsnadir_ip.so2.rtp'];

% Original R49 LAYERS RTP:
% Defaults to: satzen=0, solzen=150. emis=sea emissivity. spres = 1013.25
srcdr = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/';
srcfn = [srcdr  'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1013_400ppm.op.rtp'];

[head hatt prof patt] = rtpread(srcfn);

[nr np] = size(prof.emis);

% ----------------------------------
%     add NH3 & Zero lowest 4 layers
% ----------------------------------
nh3rtp=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
        'REGR49_400ppm_H2016_Mar2018_NH3/stdNH3_1100mb_op_400ppm.rtp'];

[hd3 ha3 pd3 pa3] = rtpread(nh3rtp);
head.glist = hd3.glist;
head.gunit = hd3.gunit;
head.ngas  = hd3.ngas;

prof.gas_11 = pd3.gas_11;
%prof.gas_11(98:101,:) = 0.0;
prof.gtotal(9,np)  = -9999;
prof.gxover(9,:)   = prof.gxover(8,:);  % shift HNO3 to gas#9
%prof.gxover(8,:)   = prof.gxover(8,:);  % make HN3 same as HNO3.
%
prof.spres = 1013.250 * ones(1, np);
prof.nlevs = 98 * ones(1, np);

% Set up variables:
scanang = [0    8.8322   17.9223   32.8244   44.8285   53.4704   59.8336 ...
        65.0428   69.3834   73.2507   76.5523   79.3679   81.7153   83.6348];

% ----------------------------------------------------------------
%   Standard Test - variable view angle, sea emissivity, night
% ----------------------------------------------------------------
h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,7);
   py.satzen(:,2)   = scanang(2);
   py.satzen(:,3)   = scanang(3);
   py.satzen(:,4)   = scanang(4);
   py.satzen(:,5)   = scanang(5);
   py.satzen(:,6)   = scanang(6);
   py.satzen(:,7)   = scanang(7);
   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
end
clear hy py;
fprintf(1,'\n')

% write out pertubation set:
outdr = '/home/chepplew/data/sarta/prod_2019/generic/';
fortp = 'r49_1013_400p_seaemis_7angs_night.rtp';

rtpwrite([outdr fortp], h2, hatt, p2, patt);

% -------------------------------------------------------------------------
% Make a version with unit emissivity. NB condition: p.rho = (1-p.emis)/pi;
% -------------------------------------------------------------------------

p3=p2;
[nemis nprof] = size(p2.emis);
p3.emis = ones(nemis,nprof);
p3.rho  = zeros(nemis,nprof);

[h4,p4] = cat_rtp(h2,p2,h2,p3);

outdr = '/home/chepplew/data/sarta/prod_2019/generic/';
fortp = 'r49_1100_98lev_400p_unitemis_seaemis_7angs_night.rtp';
rtpwrite([outdr fortp], h4, hatt, p4, patt);


