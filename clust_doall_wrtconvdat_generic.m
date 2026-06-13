addpath /home/sergio/MATLABCODE

system_slurm_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this sets tops, comment, pnums

%setting_profiles_ecm83
setting_profiles_regr49

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('do not forget to set      set_source_dir_and_rtp_files ')
disp('do not forget to set      set_source_dir_and_rtp_files ')
disp('do not forget to set      set_source_dir_and_rtp_files ')

%% kleenslurm; sbatch --array=1-49 sergio_matlab_chipCPU2024.sbatch 1
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%% JOB = 1 : 49 or 1 : 84 or 1 : 704
if length(JOB) == 0
  JOB = 1;
end

topts.myset  = 'set1';     %% do 5 set set separatelu
  [ok] = doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set2';     %% do 5 set set separatelu
  [ok] = doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set3';     %% do 5 set set separatelu
  [ok] = doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set4';     %% do 5 set set separatelu
  [ok] = doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set5';     %% do 5 set set separatelu
  [ok] = doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);

