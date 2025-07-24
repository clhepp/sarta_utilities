addpath /home/sergio/MATLABCODE

system_slurm_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this sets tops, comment, pnums
setting_profiles_ecm83
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%% JOB = 1 : 240
if length(JOB) == 0
  JOB = 1;
end

topts.myset  = 'set1';     %% do 5 set set separatelu
  [ok]=doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set2';     %% do 5 set set separatelu
  [ok]=doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set3';     %% do 5 set set separatelu
  [ok]=doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set4';     %% do 5 set set separatelu
  [ok]=doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);
topts.myset  = 'set5';     %% do 5 set set separatelu
  [ok]=doall_wrtconvdat_generic(topts, comment, pnums, [], JOB);

