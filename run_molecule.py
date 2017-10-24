
import time
import simtk.openmm as mm
from simtk.unit import *
import numpy as np
from mdtraj.reporters import NetCDFReporter
from smarty import *
import sys
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, generateTopologyFromOEMol, read_molecules
from sys import stdout
from solvationtoolkit.solvated_mixtures import *
import os

cwd = os.getcwd()
"""
# Set up mixture of ethanol and cyclohexane
# using specified numbers of each component
# solvationtoolkit uses openmoltools to build appropriately sized box for you!!!
mixture = MixtureSystem(cwd)
mixture.addComponent(label=sys.argv[1], smiles=sys.argv[2], number=int(sys.argv[3]))
#mixture.addComponent(label=sys.argv[4], smiles=sys.argv[5], number=int(sys.argv[6]))
#Generate output files for AMBER
mixture.build()
"""

num_steps = 2500000
traj_freq = 1000
data_freq = 1000
mol2files = ['monomers/'+sys.argv[1]+'.mol2']#,'monomers/'+sys.argv[4]+'.mol2']

flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
mols = []
mol = oechem.OEMol()
for mol2file in mol2files:
    ifs = oechem.oemolistream(mol2file)
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        oechem.OETriposAtomNames(mol)
        mols.append(oechem.OEGraphMol(mol))

name = sys.argv[1]+'_'+sys.argv[3]
filename = 'packmol_boxes/'+name+'.pdb'
pdb = PDBFile(filename)
forcefield = ForceField(get_data_filename('forcefield/smirnoff99Frosst.ffxml'))

params = forcefield.getParameter(smirks='[#1:1]-[#8]')
params['rmin_half']='0.01'
params['epsilon']='0.01'
forcefield.setParameter(params, smirks='[#1:1]-[#8]')

smirkseries = sys.argv[4]
eps = sys.argv[5]
rmin = sys.argv[6]
epsval = sys.argv[7]
rminval = sys.argv[8]

param = forcefield.getParameter(smirks=smirkseries)
param[eps] = epsval
param[rmin] = rminval
forcefield.setParameter(param, smirks=smirkseries)

system = forcefield.createSystem(pdb.topology,mols,nonbondedMethod=PME,nonbondedCutoff=1.2*nanometers)
integrator = LangevinIntegrator(293*kelvin, 1/picosecond, 0.002*picoseconds)
barostat = MonteCarloBarostat(1.01*bar, 293.0*kelvin, 25)
system.addForce(barostat)

#platform = mm.Platform.getPlatformByName('Reference')
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed','DeterministicForces': 'true'}

simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

netcdf_reporter = NetCDFReporter('traj_cychex_neat/'+name+'_'+smirkseries+'_'+eps+epsval+'_'+rmin+rminval+'.nc', traj_freq)
pdb_reporter = PDBReporter('traj_cychex_neat/'+name+'_'+smirkseries+'_'+eps+epsval+'_'+rmin+rminval+'.pdb', traj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(StateDataReporter(stdout, data_freq, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(StateDataReporter('StateData_cychex_neat/'+name+'_'+smirkseries+'_'+eps+epsval+'_'+rmin+rminval+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

print("Starting simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")

#Do simulation
#integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
#platform = mm.Platform.getPlatformByName('Reference')
#simulation = app.Simulation(topology, system, integrator)
#simulation.context.setPositions(positions)
#simulation.context.setVelocitiesToTemperature(temperature*kelvin)
#netcdf_reporter = NetCDFReporter('traj/'+molname+'.nc', trj_freq)
#simulation.reporters.append(netcdf_reporter)
#simulation.reporters.append(app.StateDataReporter('StateData/data_'+molname+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

#print("Starting simulation")
#start = time.clock()
#simulation.step(num_steps)
#end = time.clock()

#print("Elapsed time %.2f seconds" % (end-start))
#netcdf_reporter.close()
#print("Done!")
