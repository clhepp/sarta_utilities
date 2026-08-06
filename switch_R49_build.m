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
    
  case 'june2026_regr49_pbl'    % AIRS_OCO2_PBL new LAYER new LAYER version, CRIS
    rdpath=['/umbc/xfs3/strow/sergio_test/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
		'REGR49_400ppm_H2024_Jun2026_PBL_AIRS2834_3CrIS_IASI/'];
    rtpfile   = [rdpath 'us_std_for_pbl_breakouts400.op.rtp'];      
    dpath = [rdpath '/BREAKOUTS/'];
    iFoundBuild = +1;
    
  case 'june2026_regr49'    % AIRS 100  usual LAYER version, CRIS
    rdpath=['/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
            'REGR49_400ppm_H2024_Mar2026_AIRS2834_3CrIS_IASI/'];
    rtpfile   = [rdpath 'regr49_1100_400ppm_unitemiss.op.rtp'];
    dpath = [rdpath '/BREAKOUTS/'];
    iFoundBuild = +1;
    
  case 'aug2026_regr49_aircraft_12km'    % 12 km aircraft, AIRS
    rdpath  = ['/umbc/xfs3/strow/sergio_test/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2024_Aug2026_AIRCRAFT_12km/'];
    rtpfile = [rdpath 'regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_12km.op.rtp'];
    dpath = [rdpath '/BREAKOUTS/'];
    iFoundBuild = +1;
    
  case 'aug2026_regr49_aircraft_12km_EXPT'    % 12 km aircraft, AIRS
    rdpath  = ['/umbc/xfs3/strow/sergio_test/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2024_Aug2026_AIRCRAFT_12km/'];
    rtpfile = [rdpath 'regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_12km.op.rtp'];
    dpath = [rdpath '/BREAKOUTS/'];
    iFoundBuild = +1;
    
  case 'aug2026_regr49_aircraft_20km'    % 20 km aircraft, AIRS
    rdpath  = ['/umbc/xfs3/strow/sergio_test/git/matlabcode/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2024_Aug2026_AIRCRAFT_20km/'];
    rtpfile = [rdpath 'regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_20km.op.rtp'];
    dpath = [rdpath '/BREAKOUTS/'];
    iFoundBuild = +1;
end
