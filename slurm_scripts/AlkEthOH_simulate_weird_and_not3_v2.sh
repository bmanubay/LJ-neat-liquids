#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 12:00:00
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
python run_molecule_static_v2.py AlkEthOH_c856 &
python run_molecule_static_v2.py AlkEthOH_c854 &
python run_molecule_static_v2.py AlkEthOH_c650 &
python run_molecule_static_v2.py AlkEthOH_c598 &
python run_molecule_static_v2.py AlkEthOH_c853 &
python run_molecule_static_v2.py AlkEthOH_c841 &
python run_molecule_static_v2.py AlkEthOH_c663 &
python run_molecule_static_v2.py AlkEthOH_c1298 &
python run_molecule_static_v2.py AlkEthOH_c1163 &
python run_molecule_static_v2.py AlkEthOH_c1143 &
python run_molecule_static_v2.py AlkEthOH_c1008 &
python run_molecule_static_v2.py AlkEthOH_c901 &
python run_molecule_static_v2.py AlkEthOH_c1162 &
python run_molecule_static_v2.py AlkEthOH_c1285 &
python run_molecule_static_v2.py AlkEthOH_r1074 &
python run_molecule_static_v2.py AlkEthOH_r480 &
python run_molecule_static_v2.py AlkEthOH_r439 &
python run_molecule_static_v2.py AlkEthOH_r210 &
python run_molecule_static_v2.py AlkEthOH_r135 &
python run_molecule_static_v2.py AlkEthOH_r273 &
python run_molecule_static_v2.py AlkEthOH_r1125 &
python run_molecule_static_v2.py AlkEthOH_r858 &
python run_molecule_static_v2.py AlkEthOH_r161 &
python run_molecule_static_v2.py AlkEthOH_r86 &
python run_molecule_static_v2.py AlkEthOH_r1056 &
python run_molecule_static_v2.py AlkEthOH_r155 &
python run_molecule_static_v2.py AlkEthOH_r845 &
python run_molecule_static_v2.py AlkEthOH_r58 &
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen
