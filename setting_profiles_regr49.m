%lfcow2fow = never used anymore
lfcow2fow = [];

%%%%%%%%%%%%%%%%%%%%%%%%%
topts.csens    = 'airs_l1c';
topts.prod     = '2025';
topts.build    = 'july2025_ecm83';
topts.regset   = 'ecm83';
topts.ftc_home = '/asl/s1/sergio/alldata/ftcprod/'; %FTCHOME = run_create_new_production_directories.sh

%%%%%%%%%%%%%%%%%%%%%%%%%
topts.csens    = 'airs_l1c';
topts.prod     = '2026';
topts.build    = 'apr2026';
topts.regset   = 'r49';
topts.ftc_home = '/home/sergio/nogit/ftcprod/'; %% FTCHOME in run_create_new_production_directories.sh

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.nscang = 14; %% nowadays
                   %% 12 originally first 6 zenith for LW, bands 1,2,3, last  6 for SW bands 4,5,6,7

pnums = 1:48;  % (49 is the 49 the US STD)
pnums = 1:703; % (704 is the US STD)
pnums = 1:83;  % (none of them are US STD, is this a mistake, I should have had 84??? we will find out)

comment = 'first test of Sergio running SARTA fits';

topts.myset  = 'set1';     %% do 5 set set separately
  % 2834 chans, into 7 sets of channels with no overlaps
  % band1,2,3,4,5,6,7 some channels in MW may appear in set1 and others in set2
  %   ie they are not contiguous blocaks
  % but set1,2,3,4,5 cover the whole 2834 channels

