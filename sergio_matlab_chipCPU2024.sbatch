#!/bin/bash

# run this with   sbatch --array=1-83 sergio_matlab_chipCPU2024.sbatch 1

#SBATCH --job-name=doall_wrtconvdat   ## name
#SBATCH -N1                           ## number of nodes
#SBATCH --cpus-per-task 1             ## tasks per nod

#SBATCH --partition=2024
#SBATCH --cluster=chip-cpu
#SBATCH --mem-per-cpu=18000 
#SBATCH --account=pi_strow
#SBATCH --qos=shared
#SBATCH --time=0:59:00 

if [[ "$1" -eq "" ]]; then
  # this is individual clear/cld runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_doall_wrtconvdat_generic; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is individual clear/cld runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_doall_wrtconvdat_generic; exit"
fi

