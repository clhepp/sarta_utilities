% lfcow2fow = never used anymore
% lfcow2fow = never used anymore
% lfcow2fow = never used anymore

lfcow2fow = [];

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.satelliteORaircraft = +1;  %% this is satellite at 705 km
topts.csens    = 'airs_l1c';
topts.prod     = '2025';
topts.build    = 'july2025_ecm83';
topts.regset   = 'ecm83';
topts.ftc_home = '/asl/s1/sergio/alldata/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh
comment = 'first test of Sergio running SARTA fits';
iRegrSetType = 83;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.satelliteORaircraft = +1;  %% this is satellite at 705 km
topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'apr2026_regr49';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh
comment = 'Apr 2026 : running SARTA for H2024, with new LBLRTM12.17 and CKD4.3';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.satelliteORaircraft = +1;  %% this is satellite at 705 km
topts.csens    = 'cris_hr';
topts.prod     = '2026';
topts.build    = 'june2026_regr49_pbl';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh       -------->>> this should have been /umbc/xfs3/strow/sergio_test/nogit/?????
comment = 'June 2026 : running SARTA for H2024, LBLRTM12.17, CKD4.3, PBL SARTA CRIS';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%
topts.satelliteORaircraft = +1;  %% this is satellite at 705 km
topts.csens    = 'cris_hr';
topts.prod     = '2026';
topts.build    = 'june2026_regr49';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh
comment = 'June 2026 : running SARTA for H2024, LBLRTM12.17, CKD4.3, PBL SARTA CRIS';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.satelliteORaircraft = -1;  %% this is aircraft at 20 km but this really made no difference so can set to +1
topts.satelliteORaircraft = +1;  %% this is aircraft at 20 km but this really made no difference so can set to +1
topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'aug2026_regr49_aircraft_20km';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh
comment = 'Aug 2026 : running 20 km aircraft SARTA fits';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.satelliteORaircraft = +1;  %% this is aircraft at 20 km but this really made no difference so can set to +1
topts.satelliteORaircraft = -1;  %% this is aircraft at 12 km but this really made no difference so can set to +1
topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'aug2026_regr49_aircraft_12km';     
topts.build    = 'aug2026_regr49_aircraft_12km_EXPT';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh
comment = 'Aug 2026 : running 12 km aircraft SARTA fits';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_nscanang_pnums

topts.iWriteMat = -1;  %% default, just write f77 binary file, default
topts.iWriteMat = +1;  %% debug,   also write mat file  --- script reader/plotter = plot_l2s_sets.m
