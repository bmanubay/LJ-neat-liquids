#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node 25
#SBATCH -t 15:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brma3379@colorado.edu

# echo commands to stdout 
set -x

# move to working directory
cd /pylon5/ct4s8bp/brma3379

# copy input file to working directory
cp /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/run_molecule.py .
cp -rf /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/traj_EtOH_cychex_mix .
cp -rf /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/StateData_EtOH_cychex_mix .

#commands
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.20 3.2 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.225 3.2 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.25 3.2 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.275 3.2 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.30 3.2 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.20 3.6 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.225 3.6 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.25 3.6 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.275 3.6 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.30 3.6 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.20 4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.225 4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.25 4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.275 4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.30 4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.20 4.4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.225 4.4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.25 4.4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.275 4.4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.30 4.4 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.20 4.8 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.225 4.8 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.25 4.8 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.275 4.8 &
python run_molecule.py ethanol CCO 75 cyclohexane C1CCCCC1 175 [*:1] epsilon rmin_half 0.30 4.8 &
wait

# copy directory to persistent space
cp -rf traj_EtOH_cychex_mix /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/
cp -rf StateData_EtOH_cychex_mix /pylon2/ct4s8bp/brma3379/git_repos/open-forcefield-tools/single-molecule-property-generation/
