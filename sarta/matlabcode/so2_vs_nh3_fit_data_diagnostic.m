% SO2 vs NH3 SARTA fitting

%{
% The following are the kCARTA b/o gas weights used in the L2S O/Ds.
prefix = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018/';
pre_p1 = 'so2bandF/convolved_kcarta_F_';
caVers:
    include_param: 'v1.20 17-11-23 template_kcartaV120_400_H2016.param                              '
              ckd: -1
          comment: '1   -1    1.0    3 1 0.0  103 0.0 3 0.0 !all gases except H2O and O3 have weight 1.0 (F)

pre_p2 = 'so2bandS/convolved_kcarta_so2bandS_'
    include_param: 'v1.20 17-11-23 template_kcartaV120_400_H2016.param                              '
              ckd: -1
          comment: '1   -1    0.0    1 9 1.0 !all gases have weight 0.0 except SO2=1 and water DOCOMMENT ozone=0

prefix = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018_NH3/';
pre_p1 = 'n2ohno3bandF/convolved_kcarta_F_';
    include_param: 'v1.20 18-03-26 template_kcartaV120_400_H2016_NLTEH2016.param                    '
              ckd: -1
          comment: '1   -1    1.0    3 1 0.0  103 0.0 3 0.0 !all gases except H2O and O3 have weight 1.0 (F)
	  
pre_p2 = 'nh3bandNH3/convolved_kcarta_nh3bandNH3_';
    include_param: 'v1.20 18-03-26 template_kcartaV120_400_H2016_NLTEH2016.param                    '
              ckd: -1
          comment: '1   -1    0.0    1 11 1.0 !all gases have weight 0.0 except NH3=1

%}	  
%
%
% 

addpath /asl/matlib/h4tools
	  
cd /home/chepplew/projects/sarta/prod_2018/cris_hr/

% ---------------------------
% Check refprof:
% ---------------------------
fnref = '/home/chepplew/data/sarta/prod_2018/cris_hr_mar18/Coef/refprof_nh3';
inFH=fopen(fnref,'r');
clear hdr;
for i=1:9 hdr{i} = fgetl(inFH); end
A = fscanf(inFH, '%e  %e %e %e %e %e %e %e %e %e %e %e %e %e',[14,Inf]);
% cols: 2=ALT, 3=THICK, 4=PRESS, 10=CH4, 11=SO2, 12=HNO3, 13=N2O, 14=NH3
figure(1);clf;plot(A(14,:),A(1,:),'.-', A(11,:),A(1,:),'.-'); grid on;
  xlabel('Gas amnt. kmol.cm^{-2}');ylabel('level');title('reference profile Gas Amount')
  legend('NH3','SO2');
  %saveas(gcf,'./figs/refprof_nh3_so2_amnt.png','png')

% ----------------------------------
Get 49th profile from regression set
% ----------------------------------
rtp49 = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP',...
         '/stdNH3_1100mb_op_400ppm.rtp'];
[hd ha pd pa] = rtpread(rtp49);



dfn='/home/chepplew/data/sarta/prod_2018/p400/SO2/cris_hrg4_so2_data_long.mat';
so2 = load(dfn);
so2.dfn = dfn;


dfn='/home/chepplew/data/sarta/prod_2018/p400/NH3/cris_hrg4_nh3_data_long_r49.mat';
nh3 = load(dfn);
nh3.dfn = dfn;

% Choose first profile, nadir, near SFC layer
nf1=figure(1);plot(nh3.fchan, squeeze(nh3.tauz(:,98,1)),'-',...
   nh3.fchan, squeeze(nh3.tauz11(:,98,1)),'-');

nf2=figure(2);plot(so2.fchan, squeeze(so2.tauz(:,98,1)),'-',...
   so2.fchan, squeeze(so2.tauz9(:,98,1)),'-')

% NH3 lines look too strong so reduce by x100
nh3.od11 = -log(nh3.tauz11).*0.01;
