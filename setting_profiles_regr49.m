% lfcow2fow = never used anymore
% lfcow2fow = never used anymore
% lfcow2fow = never used anymore

lfcow2fow = [];

%%%%%%%%%%%%%%%%%%%%%%%%%
topts.csens    = 'airs_l1c';
topts.prod     = '2025';
topts.build    = 'july2025_ecm83';
topts.regset   = 'ecm83';
topts.ftc_home = '/asl/s1/sergio/alldata/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh
comment = 'first test of Sergio running SARTA fits';
iRegrSetType = 83;

%%%%%%%%%%%%%%%%%%%%%%%%%
topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'apr2026_regr49';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh
comment = 'Apr 2026 : running SARTA for H2024, with new LBLRTM12.17 and CKD4.3';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.csens    = 'cris_hr';
topts.prod     = '2026';
topts.build    = 'june2026_regr49_pbl';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh       -------->>> this should have been /umbc/xfs3/strow/sergio_test/nogit/?????
comment = 'June 2026 : running SARTA for H2024, LBLRTM12.17, CKD4.3, PBL SARTA CRIS';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.csens    = 'cris_hr';
topts.prod     = '2026';
topts.build    = 'june2026_regr49';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh
comment = 'June 2026 : running SARTA for H2024, LBLRTM12.17, CKD4.3, PBL SARTA CRIS';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'aug2026_regr49_aircraft_12km';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh
comment = 'Aug 2026 : running 12 km aircraft SARTA fits';
iRegrSetType = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topts.nscang = 14; %% nowadays
                   %% 12 originally first 6 zenith for LW, bands 1,2,3, last  6 for SW bands 4,5,6,7

if iRegrSetType == 703
  pnums = 1:703; % (704 is the US STD)
elseif iRegrSetType == 83  
  pnums = 1:83;  % (none of them are US STD, is this a mistake, I should have had 84??? we will find out)
elseif iRegrSetType == 49  
  pnums = 1:48;  % (49 is the 49 the US STD)
else
  error('unknown iRegrSetType')
end

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.myset  = 'set1';     %% do 5 set set separately
  % 2834 chans, into 7 sets of channels with no overlaps
  % band1,2,3,4,5,6,7 some channels in MW may appear in set1 and others in set2
  %   ie they are not contiguous blocaks
  % but set1,2,3,4,5 cover the whole 2834 channels

