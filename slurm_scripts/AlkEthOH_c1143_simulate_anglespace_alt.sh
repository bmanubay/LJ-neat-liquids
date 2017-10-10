#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 24
#SBATCH -t 08:00:00
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
python run_molecule_var.py AlkEthOH_c1143 70.0 98.5 &
python run_molecule_var.py AlkEthOH_c1143 100.0 98.5 &
python run_molecule_var.py AlkEthOH_c1143 130.0 98.5 &
python run_molecule_var.py AlkEthOH_c1143 80.0 102.2 &
python run_molecule_var.py AlkEthOH_c1143 100.0 102.2 &
python run_molecule_var.py AlkEthOH_c1143 120.0 102.2 &
python run_molecule_var.py AlkEthOH_c1143 90.0 105.9 &
python run_molecule_var.py AlkEthOH_c1143 100.0 105.9 &
python run_molecule_var.py AlkEthOH_c1143 110.0 105.9 &
python run_molecule_var.py AlkEthOH_c1143 70.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 80.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 90.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 110.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 120.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 130.0 109.5 &
python run_molecule_var.py AlkEthOH_c1143 90.0 113.2 &
python run_molecule_var.py AlkEthOH_c1143 100.0 113.2 &
python run_molecule_var.py AlkEthOH_c1143 110.0 113.2 &
python run_molecule_var.py AlkEthOH_c1143 80.0 116.9 &
python run_molecule_var.py AlkEthOH_c1143 100.0 116.9 &
python run_molecule_var.py AlkEthOH_c1143 120.0 116.9 &
python run_molecule_var.py AlkEthOH_c1143 70.0 120.5 &
python run_molecule_var.py AlkEthOH_c1143 100.0 120.5 &
python run_molecule_var.py AlkEthOH_c1143 130.0 120.5 &
wait

# copy directory to persistent space
cp *.nc /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/trajs_for_simple_surr_gen_ang
