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
python run_molecule_static_v2.py AlkEthOH_c651 &
python run_molecule_static_v2.py AlkEthOH_c635 &
python run_molecule_static_v2.py AlkEthOH_c674 &
python run_molecule_static_v2.py AlkEthOH_c44 &
python run_molecule_static_v2.py AlkEthOH_c613 &
python run_molecule_static_v2.py AlkEthOH_c102 &
python run_molecule_static_v2.py AlkEthOH_c661 &
python run_molecule_static_v2.py AlkEthOH_c1137 &
python run_molecule_static_v2.py AlkEthOH_c1299 &
python run_molecule_static_v2.py AlkEthOH_c1144 &
python run_molecule_static_v2.py AlkEthOH_c1266 &
python run_molecule_static_v2.py AlkEthOH_c1119 &
python run_molecule_static_v2.py AlkEthOH_c1296 &
python run_molecule_static_v2.py AlkEthOH_c1282 &
python run_molecule_static_v2.py AlkEthOH_r552 &
python run_molecule_static_v2.py AlkEthOH_r351 &
python run_molecule_static_v2.py AlkEthOH_r434 &
python run_molecule_static_v2.py AlkEthOH_r485 &
python run_molecule_static_v2.py AlkEthOH_r231 &
python run_molecule_static_v2.py AlkEthOH_r78 &
python run_molecule_static_v2.py AlkEthOH_r1085 &
python run_molecule_static_v2.py AlkEthOH_r100 &
python run_molecule_static_v2.py AlkEthOH_r151 &
python run_molecule_static_v2.py AlkEthOH_r871 &
python run_molecule_static_v2.py AlkEthOH_r60 &
python run_molecule_static_v2.py AlkEthOH_r150 &
python run_molecule_static_v2.py AlkEthOH_r1142 &
python run_molecule_static_v2.py AlkEthOH_r93 &
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen
