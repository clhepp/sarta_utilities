%% this should not really need fixing unless you use eg 25000 or 38000 profile set ...
%% common code used by setting_profiles_regr49.m
%%                     setting_profiles_ecm83.m

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
  error('set_nscanang_pnums.m : unknown iRegrSetType')
end

%%%%%%%%%%%%%%%%%%%%%%%%%

topts.myset  = 'set1';     %% do 5 set set separately
  % 2834 chans, into 7 sets of channels with no overlaps
  % band1,2,3,4,5,6,7 some channels in MW may appear in set1 and others in set2
  %   ie they are not contiguous blocaks
  % but set1,2,3,4,5 cover the whole 2834 channels

