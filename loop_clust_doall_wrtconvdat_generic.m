JMAX = 49;   %% number of profiles, including US STD

JMAX = JMAX - 1;

for JOBB = 1 : JMAX
  clear JOB
  clust_doall_wrtconvdat_generic
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% should be JMAX x 5 = 48x5 = 240 files for 49 regr
watch "ls -lt /home/sergio/git/ftc_dev/outd/../*/*.dat | wc -l"
%}
