#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 25:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

# echo commands to stdout 
set -x

# move to working directory
cd /pylon5/ct4s8bp/brma3379

# copy input file to working directory
#cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/run_molecule_static.py .
#cp -rf /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/Mol2_files .

#commands
python run_molecule_static.py AlkEthOH_c603 &
python run_molecule_static.py AlkEthOH_c844 &
python run_molecule_static.py AlkEthOH_c88 &
python run_molecule_static.py AlkEthOH_c653 &
python run_molecule_static.py AlkEthOH_c600 &
python run_molecule_static.py AlkEthOH_c845 &
python run_molecule_static.py AlkEthOH_c699 &
python run_molecule_static.py AlkEthOH_c1160 &
python run_molecule_static.py AlkEthOH_c24 &
python run_molecule_static.py AlkEthOH_c1295 &
python run_molecule_static.py AlkEthOH_c1283 &
python run_molecule_static.py AlkEthOH_c1135 &
python run_molecule_static.py AlkEthOH_c258 &
python run_molecule_static.py AlkEthOH_c488 &
python run_molecule_static.py AlkEthOH_r368 &
python run_molecule_static.py AlkEthOH_r639 &
python run_molecule_static.py AlkEthOH_r397 &
python run_molecule_static.py AlkEthOH_r392 &
python run_molecule_static.py AlkEthOH_r417 &
python run_molecule_static.py AlkEthOH_r268 &
python run_molecule_static.py AlkEthOH_r363 &
python run_molecule_static.py AlkEthOH_r13 &
python run_molecule_static.py AlkEthOH_r3 &
python run_molecule_static.py AlkEthOH_r1136 &
python run_molecule_static.py AlkEthOH_r1012 &
python run_molecule_static.py AlkEthOH_r181 &
python run_molecule_static.py AlkEthOH_r166 &
python run_molecule_static.py AlkEthOH_r68 &
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen
