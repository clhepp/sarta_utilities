iFoundBuild = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(regset,'R49'))
  fprintf(1,' looking for regset = R49  build = %s \n',build)
  switch_R49_build
  
  pnums     = [1:48]';
  comment   = [csens ' r49 400ppm H2020 ftc.14a ' build];
%  outd_pref = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' build '/'];
  outd_pref = [FTC_HOME prod_run '/' lower(csens) '/' build '/'];

  if iFoundBuild < 0
    error('set_source_dir_and_rtp_files.m : strcmp(regset,R49) and iFoundBuild < 0')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  if iFoundBuild < 0
    error('set_source_dir_and_rtp_files.m : strcmp(regset,SAF704) and iFoundBuild < 0')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  if iFoundBuild < 0
    error('set_source_dir_and_rtp_files.m : strcmp(regset,ECM83) and iFoundBuild < 0')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iFoundBuild == -1
  fprintf(1,'regset = %s build = %s combo \n',regset,build)
  error('DID NOT FIND this regset/build combo ::: please check doall_wrtconvdat_generic.m')
end
