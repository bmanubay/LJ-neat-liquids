#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node 16
#SBATCH -t 30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

source activate my_root

# echo commands to stdout 
set -x

# move to working directory
cd /pylon5/ct4s8bp/brma3379

# copy input file to working directory
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/run_molecule_static.py .
#cp -rf /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/Mol2_files .

#commands
python run_molecule_static.py AlkEthOH_c636 &
python run_molecule_static.py AlkEthOH_c868 &
python run_molecule_static.py AlkEthOH_c821 &
python run_molecule_static.py AlkEthOH_c820 &
python run_molecule_static.py AlkEthOH_c1284 &
python run_molecule_static.py AlkEthOH_c38 &
python run_molecule_static.py AlkEthOH_c0 &
python run_molecule_static.py AlkEthOH_c1302 &
python run_molecule_static.py AlkEthOH_r574 &
python run_molecule_static.py AlkEthOH_r82 &
python run_molecule_static.py AlkEthOH_r283 &
python run_molecule_static.py AlkEthOH_r1138 &
python run_molecule_static.py AlkEthOH_r1152 &
python run_molecule_static.py AlkEthOH_r182 &
python run_molecule_static.py AlkEthOH_r51 &
python run_molecule_static.py AlkEthOH_r131 &
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen
