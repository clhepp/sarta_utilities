iFoundBuild = -1;
if (strcmp(regset,'R49'))
  fprintf(1,' looking for regset = R49  build = %s \n',build)
  switch build
    case 'jun2016'
      kpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
               'REGR49_400ppm_H2012_June2016/'];
      iFoundBuild = +1;      
    case 'mar2018'
     %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/';  
     dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
              'REGR49_400ppm_H2016_Mar2018/'];
      iFoundBuild = +1;           
    case 'sep2018'  
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
               'REGR49_400ppm_H2016_Sept2018_AIRS2645/'];
      iFoundBuild = +1;            
    case 'dec2018' 
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
               'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];
      iFoundBuild = +1;            
    case 'feb2020'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
               'REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/'];
      iFoundBuild = +1;            
    case 'may2021'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
               'REGR49_400ppm_H2016_May2021_AIRS2834_3CrIS_IASI/'];
      iFoundBuild = +1;            
    case 'jul2022'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [dpath 'regr49_1100_400ppm_unitemiss.op.rtp'];
      iFoundBuild = +1;            
    case 'jan2025a'    % AIRS_OCO2_PBL new LAYER version
      dpath=['/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2020_Jan2025_PBL_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [dpath 'regr49_pbl.op.rtp'];
      iFoundBuild = +1;            
    case 'apr2026'    % AIRS 100  usual LAYER version
      rdpath=['/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
		'REGR49_400ppm_H2024_Mar2026_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [rdpath 'regr49_1100_400ppm_unitemiss.op.rtp'];      
      dpath = [rdpath '/BREAKOUTS/'];
      iFoundBuild = +1;            
    case 'june2026_regr49_pbl'    % AIRS_OCO2_PBL new LAYER new LAYER version
      rdpath=['/umbc/xfs3/strow/sergio_test/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
		'REGR49_400ppm_H2024_Jun2026_PBL_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [rdpath 'us_std_for_pbl_breakouts400.op.rtp'];      
      dpath = [rdpath '/BREAKOUTS/'];
      iFoundBuild = +1;            
  end
  pnums     = [1:48]';
  comment   = [csens ' r49 400ppm H2020 ftc.14a ' build];
%  outd_pref = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' build '/'];
  outd_pref = [FTC_HOME prod_run '/' lower(csens) '/' build '/'];
end

if(strcmp(regset,'SAF704'))
  %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/';
  dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
         'SAF704_400ppm_H2016_Dec2018_AIRS2834/'];
  rtpfile = [dpath 'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
  rtpfile = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/' ...
             'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
  pnums = [1:703]';
  comment = [csens ' SAF704 400ppm CO2 H2016'];
%  outd_pref = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' build '/'];
  outd_pref = [FTC_HOME prod_run '/' lower(csens) '/' build '/'];
  iFoundBuild = +1;        
end

if(strcmp(regset,'ECM83'))
  %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/';
  dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
           'REGR49_400ppm_H2020_July2025_ECMWF83Profiles_AIRS2834_3CrIS_IASI/'];
  rtpfile = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/ECMWF_83P_91L_Sept2015/' ...
               'ecmwf_co2_400ppm_1100mb.op.rtp'];
  pnums = [1:83]';
  comment = [csens ' ECM83 400ppm CO2 H2020'];
  outd_pref = [FTC_HOME prod_run '/' lower(csens) '/' build '/'];
  iFoundBuild = +1;          
end

if iFoundBuild == -1
  fprintf(1,'DID NOT FINDver found the regset = %s build = %s combo \n',regset,build)
  error('please check doall_wrtconvdat_generic.m')
end
