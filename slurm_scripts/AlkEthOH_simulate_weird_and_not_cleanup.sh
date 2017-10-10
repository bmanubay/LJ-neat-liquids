#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node 25
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

# echo commands to stdout 
set -x

# move to working directory
cd /pylon5/ct4s8bp/brma3379

# copy input file to working directory
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/run_molecule_var.py .
#cp -rf /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/Mol2_files .

#commands
python run_molecule_var.py AlkEthOH_c868 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c1008 980.0 0.8 &
python run_molecule_var.py AlkEthOH_r283 980.0 0.8 &
python run_molecule_var.py AlkEthOH_r845 980.0 0.8 &
python run_molecule_var.py AlkEthOH_r51 980.0 0.8 &
python run_molecule_var.py AlkEthOH_r480 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c1284 680.0 1.09 &
python run_molecule_var.py AlkEthOH_r1152 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c820 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c0 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c1302 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c821 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c636 980.0 0.8 &
python run_molecule_var.py AlkEthOH_c868 680.0 1.4 & 
python run_molecule_var.py AlkEthOH_r82 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_r166 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_r1138 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_r161 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_c1285 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_c845 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_r1056 680.0 1.09 & 
python run_molecule_var.py AlkEthOH_r182 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_r131 680.0 1.09 & 
python run_molecule_var.py AlkEthOH_c853 980.0 0.8 & 
python run_molecule_var.py AlkEthOH_c38 980.0 0.8 & 
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen
