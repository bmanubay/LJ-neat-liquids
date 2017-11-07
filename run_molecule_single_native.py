#!/bin/env python

import time
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *
import numpy as np
from mdtraj.reporters import NetCDFReporter
# Import the SMIRNOFF forcefield engine and some useful tools
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye import oechem
import sys
import numpy as np

molname = [sys.argv[1]]
mol_filename = ['monomers/'+m+'.mol2' for m in molname]
time_step = 0.8 #Femtoseconds
temperature = 293.15 #kelvin
friction = 1 # per picosecond
num_steps = 7500000
trj_freq = 1000 #steps
data_freq = 1000 #steps

# Load OEMol
for ind,j in enumerate(mol_filename):
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(j)
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    oechem.OETriposAtomNames(mol)

    # Get positions
    coordinates = mol.GetCoords()
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float64)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = Quantity(positions, angstroms)
    
    
    # Load forcefield
    forcefield = ForceField(get_data_filename('forcefield/smirnoff99Frosst.ffxml'))

    # Define system
    topology = generateTopologyFromOEMol(mol)
    params = forcefield.getParameter(smirks='[#1:1]-[#8]')
    params['rmin_half']='0.01'
    params['epsilon']='0.01'
    forcefield.setParameter(params, smirks='[#1:1]-[#8]')
    system = forcefield.createSystem(topology, [mol])
    """
    paramlist1 = np.array([100.,140.,60.,140.,60.])#np.arange(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
    paramlist2 = np.array([109.5,76.5,76.5,142.5,142.5])
    #j = sys.argv[5]
    
    smirkseries = '[*:1]~[#6X4:2]-[*:3]'#sys.argv[6]#'[#6X4:1]-[#1:2]'
    
    paramtype1 = 'k'#sys.argv[7]#'length'
    paramtype2 = 'angle'#sys.argv[8]


    param = forcefield.getParameter(smirks=smirkseries)
    for i,j in zip(paramlist1,paramlist2):
        param[paramtype1] = str(i)
        param[paramtype2] = str(j)
        forcefield.setParameter(param, smirks=smirkseries)
        system = forcefield.createSystem(topology, [mol])

    """
    #Do simulation
    integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed','DeterministicForces': 'true'}
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature*kelvin)
    netcdf_reporter = NetCDFReporter(molname[ind]+'_vacuum.nc', trj_freq)#NetCDFReporter('traj4ns_c1143/'+molname[ind]+'_'+smirkseries+'_'+paramtype1+str(i)+'_'+paramtype2+str(j)+'.nc', trj_freq)
    simulation.reporters.append(netcdf_reporter)
        #simulation.reporters.append(app.StateDataReporter('StateData4ns_c1143/data_'+molname[ind]+'_'+smirkseries+'_'+paramtype1+str(i)+'_'+paramtype2+str(j)+'.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
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
