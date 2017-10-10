#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node 3 
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

# echo commands to stdout 
set -x

# move to working directory
cd /pylon5/ct4s8bp/brma3379

# copy input file to working directory
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/graph_bl_test_*.py .
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/AlkEthOH_test_C-H_len_evidence.csv .
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/AlkEthOH_test_C-C-C_ang_evidence.csv .

#commands
python graph_bl_test_simp.py &
python graph_bl_test_1.py &
python graph_bl_test_1int.py &
wait

# copy directory to persistent space
cp *.csv /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation
