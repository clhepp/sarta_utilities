# Collect atmospheric profiles for extended nonLTE application
  Data sources are: SAF704, REGR49, SABER, AFGL, WACCM-X 
   Matlab code scripts are:
      concat_saf_and_regr49.m
     splice_waccm_to_saf.m
     load_regr49_to_147levs.m
     subset_saf.m
     load_regr49_add_afgl_saber.m
     subset_saf_append_r49_afgl_saber.m
     splice_afgl_to_saf_r49.m
     append_r49_to_saf.m
     get_afgl_atmos_type.m
## WorkFlow
1. subset_saf  -> splice_afgl_to_saf   ->
               -> splice waccm_to_saf  ->  concat_saf_and_regr49

2. load_regr49_to_147levs.m  -> splice_afgl_to_regr49   ->
                             -> splice_waccm_to_regr49  -> concat...

# Output:
  saved to: /home/chepplew/data/sarta/SAF/saf_sub300_reg49_plv147_400ppm

# Notes:
Original SAF: 
wget https://nwp-saf.eumetsat.int/site/download/profile_datasets/profiles137.tar.bz2
at ASL: /asl/rta/ECMWF_RTTOV_91_25000_and_704_Profiles/ECMWF_SAF_137ProfilesStuff/ECMWF_SAF_137Profiles/
or: /home/chepplew/data/SAF/SAF137


