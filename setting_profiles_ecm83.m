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
topts.ftc_home = '/asl/s1/sergio/alldata/ftcprod/';
  %FTCHOME = run_create_new_production_directories.sh
comment = 'first test of Sergio running SARTA fits';  
iRegrSetType = 83;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_nscanang_pnums

topts.iWriteMat = +1;  %% debug,   also write mat file
topts.iWriteMat = -1;  %% default, just write f77 binary file, default
