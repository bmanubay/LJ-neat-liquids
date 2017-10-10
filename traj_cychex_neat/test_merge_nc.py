import xray as xr
from sys import exit
from pdb import set_trace
import netCDF4 as nc
import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import pandas as pd
import pymbar as mb
from pymbar import timeseries
from collections import OrderedDict
from smarty import *
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, generateTopologyFromOEMol, read_molecules
import mdtraj as md

#-------------------------------------------------
def read_traj(ncfiles,indkeep=0):
    """
    Take multiple .nc files and read in coordinates in order to re-valuate energies based on parameter changes

    Parameters
    -----------
    ncfiles - a list of trajectories in netcdf format

    Returns
    ----------
    data - all of the data contained in the netcdf file
    xyzn - the coordinates from the netcdf in angstroms
    """

    data = nc.Dataset(ncfiles)
    xyz = data.variables['coordinates']
    lens = data.variables['cell_lengths']
    angs = data.variables['cell_angles']
    angsn = np.float64(angs[indkeep:-2])#Quantity(angs[indkeep:-2], degrees)
    lensn = np.float64(lens[indkeep:-2])#Quantity(lens[indkeep:-2], angstroms)
    xyzn = np.float64(xyz[indkeep:-2])#Quantity(xyz[indkeep:-2].data, angstroms)
    
    return data, xyzn, lensn, angsn
#------------------------------------------------------------------
def get_energy(system, positions, vecs):
    """
    Return the potential energy.
    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    vecs : simtk.unit.Quantity of dimension 3 with unit of length
        Box vectors to use 
    Returns
    ---------
    energy
    """
   
    integrator = mm.VerletIntegrator(1.0 * femtoseconds)
    context = mm.Context(system, integrator)
    context.setPositions(positions*0.1)
    vecs = [0.1*a for a in vecs]
    vecs = np.array(vecs)
    print vecs
    context.setPeriodicBoxVectors(*vecs)
    
    state = context.getState(getEnergy=True)
    
    energy = state.getPotentialEnergy() 
    return energy

#---------------------------------------------------
def new_param_energy(coords, params, topology, vecs, P=1.01, T=300.):
    """
    Return potential energies associated with specified parameter perturbations.
    Parameters
    ----------
    coords: coordinates from the simulation(s) ran on the given molecule
    params:  arbitrary length dictionary of changes in parameter across arbitrary number of states. Highest level key is the molecule AlkEthOH_ID,
             second level of keys are the new state, the values of each of these subkeys are a arbitrary length list of length 3 lists where the
             length 3 lists contain information on a parameter to change in the form: [SMIRKS, parameter type, parameter value]. I.e. :

             params = {'AlkEthOH_c1143':{'State 1':[['[6X4:1]-[#1:2]','k','620'],['[6X4:1]-[#6X4:2]','length','1.53'],...],'State 2':[...],...}}
    P: Pressure of the system. By default set to 1.01 bar.
    T: Temperature of the system. By default set to 300 K.

    Returns
    -------
    E_kn: a kxN matrix of the dimensional energies associated with the forcfield parameters used as input
    u_kn: a kxN matrix of the dimensionless energies associated with the forcfield parameters used as input
    """

    #-------------------
    # CONSTANTS
    #-------------------
    kB = 0.0083145  #Boltzmann constant (Gas constant) in kcal/(mol*K)
    beta = 1/(kB*T)

    #-------------------
    # PARAMETERS
    #-------------------
    params = params

    # Determine number of states we wish to estimate potential energies for
    mol2files = []
    for i in params:
        mol2files.append('monomers/'+i.rsplit(' ',1)[0]+'.mol2')
    #    mol2files.append('monomers/'+i.rsplit(' ',1)[1]+'.mol2')
    #print mol2files
    #mol = 'Mol2_files/'+mols[0]+'.mol2'
    #K = len(params[mols[0]].keys())


    #if np.shape(params) != np.shape(N_k): raise "K_k and N_k must have same dimensions"


    # Determine max number of samples to be drawn from any state

    #mol2files = ['monomers/'+sys.argv[1]+'.mol2','monomers/'+sys.argv[4]+'.mol2']

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
    K = len(params['cyclohexane'].keys())
    
    # Load forcefield file
    ffxml = get_data_filename('forcefield/smirnoff99Frosst.ffxml')
    ff = ForceField(ffxml)

    # Generate a topology
    #from smarty.forcefield import generateTopologyFromOEMol
    top = topology#generateTopologyFromOEMol(mol)
     
    #-----------------
    # MAIN
    #-----------------

    # Calculate energies

    E_kn = np.zeros([K,len(coords)],np.float64)
    u_kn = np.zeros([K,len(coords)],np.float64)
    for i,j in enumerate(params):
        AlkEthOH_id = j
        for k,l in enumerate(params[AlkEthOH_id]):
            for m,n in enumerate(params[AlkEthOH_id][l]):
                newparams = ff.getParameter(smirks=n[0])
                newparams[n[1]]=n[2]
                ff.setParameter(newparams,smirks=n[0])
                system = ff.createSystem(top,mols,nonbondedMethod=PME,nonbondedCutoff=12.*angstroms)
                barostat = MonteCarloBarostat(P*bar, T*kelvin, 25)
                system.addForce(barostat)
            
                #print system
                
            for o,p in enumerate(coords):
                e = get_energy(system,p,vecs[o])
                print o,e
                E_kn[k,o] = e._value
                u_kn[k,o] = e._value*beta
                

    return E_kn,u_kn
#----------------------------------------------------
def subsampletimeseries(timeser,xyzn,N_k):
    """
    Return a subsampled timeseries based on statistical inefficiency calculations.
    Parameters
    ----------
    timeser: the timeseries to be subsampled
    xyzn: the coordinates associated with each frame of the timeseries to be subsampled
    N_k: original # of samples in each timeseries

    Returns
    ---------
    N_k_sub: new number of samples per timeseries
    ts_sub: the subsampled timeseries
    xyz_sub: the subsampled configuration series
    """
    # Make a copy of the timeseries and make sure is numpy array of floats
    ts = timeser
    xyz = xyzn

    # initialize array of statistical inefficiencies
    g = np.zeros(len(ts),np.float64)


    for i,t in enumerate(ts):
        if np.count_nonzero(t)==0:
            g[i] = np.float(1.)
            print "WARNING FLAG"
        else:
            g[i] = timeseries.statisticalInefficiency(t)

    N_k_sub = np.array([len(timeseries.subsampleCorrelatedData(t,g=b)) for t, b in zip(ts,g)])
    ind = [timeseries.subsampleCorrelatedData(t,g=b) for t,b in zip(ts,g)]

    if (N_k_sub == N_k).all():
        ts_sub = ts
        xyz_sub = xyz
        print "No sub-sampling occurred"
    else:
        print "Sub-sampling..."
        ts_sub = np.array([t[timeseries.subsampleCorrelatedData(t,g=b)] for t,b in zip(ts,g)])
        #for c in xyz:
        #    xyz_sub = [c[timeseries.subsampleCorrelatedData(t,g=b)] for t,b in zip(ts,g)]
        for i,j in enumerate(xyz):
            xyz_sub = [j[ii] for ii in ind[i]]

    return ts_sub, N_k_sub, xyz_sub, ind
#-------------------------------------------------------------------------

kB = 0.0083145
T= 293.15 
"""
#files = nc.glob('*_cont.nc')
#files = files[0]
files = ['cyclohexane_250_[#6X4:1]_epsilon0.1094_rmin_half1.9080_cont.nc']#,'cyclohexane_250_[#6X4:1]_epsilon0.05_rmin_half1.9080_cont.nc']
#print files[0]
file_strings = [i.rsplit('_',1)[0] for i in files]
file_str = set(file_strings)

file_tups_traj = [[i+'.nc',i+'_cont.nc'] for i in file_str]

#file_tups_traj = file_tups_traj[0]
file_tups_sd = [[i+'.csv',i+'_cont.csv'] for i in file_str]
#file_tups_sd = file_tups_sd[0]
params = [i.rsplit('.',1)[0].rsplit('_') for i in files]
params = [[i[3][7:],i[5][4:]] for i in params]
MMcyc = 0.08416 #kg/mol
#with md.NetCDFTrajectoryFile(files[0]) as f:
#    xyz, cell_lengths, cell_angles = f.read()
#filepdb = 'cyclohexane_250_[#6X4:1]_epsilon0.15_rmin_half1.9080.pdb'
#filename = 'cyclohexane_250.pdb'
#pdb = PDBFile(filename)
#top = pdb.topology
#print pdb, top
#a = md.load_netcdf(files[0],filename)
#a = md.formats.NetCDFTrajectoryFile('cyclohexane_250_[#6X4:1]_epsilon0.15_rmin_half1.9080.nc')#[0])
#dat = a.read()
#a = md.load_pdb(filepdb)
#set_trace()
states_traj = [[] for i in file_tups_traj] 
states_sd = [[] for i in file_tups_sd]
xyz_orig = [[] for i in file_tups_traj]
vol_orig = [[] for i in file_tups_traj]
lens_orig = [[] for i in file_tups_traj]
angs_orig = [[] for i in file_tups_traj]
ener_orig = [[] for i in file_tups_sd]
dens_orig = [[] for i in file_tups_sd]
vecs_orig = [[] for i in file_tups_sd]
burnin = 1800
print 'Analyzing Cyclohexane neat liquid trajectories'
for j,i in enumerate(file_tups_traj):
    #print 'Analyzing trajectory %s of %s'%(j+1,len(file_tups_traj)+1)
    for ii in i:    
        try:
            data, xyz, lens, angs = read_traj(ii,burnin)             
            #set_trace()
        except IndexError:
            print "The trajectory had fewer than %s frames" %(burnin)
            continue 
            #print 'Well, that one is fucked \(`_`\)'
        for m,n in zip(lens,angs):  
            vecs = md.utils.lengths_and_angles_to_box_vectors(float(m[0]),float(m[1]),float(m[2]),float(n[0]),float(n[1]),float(n[2]))
            vecs_orig[j].append(vecs)
        #print vecs_orig
        for pos in xyz:
            xyz_orig[j].append(pos)
        for l in lens:
            vol_orig[j].append(np.prod(l)) #A**3
            lens_orig[j].append(l)
        for ang in angs:
            angs_orig[j].append(ang)
    states_traj[j].append(i[0].rsplit('.',1)[0])
"""
"""
for j,i in enumerate(file_tups_sd):
    try:
        datasets = [pd.read_csv(ii,sep=',')[250:] for ii in i]
        merged = pd.concat(datasets)
    except IndexError:
        print "The state data record had fewer than 250 frames"
    for e in merged["Potential Energy (kJ/mole)"]:
        ener_orig[j].append(e)
    for dens in merged["Density (g/mL)"]:
        dens_orig[j].append(dens)
    states_sd[j].append(i[0].rsplit('.',1)[0])

# Compute enthalpies H = U + pV
p = 101000. #Pa

H = [[] for i in file_tups_sd]

for j,i in enumerate(vol_orig):
    for jj,ii in enumerate(i):
        try: 
            H[j].append(ener_orig[j][jj]*1.e3 + p*ii*1.e-30)
        except IndexError:
            print "the indices don't match for some reason"


print states_traj
print len(states_traj)
print states_sd
print len(states_sd)
print [len(i) for i in xyz_orig]
print len(xyz_orig)
print [len(i) for i in vol_orig]
print len(vol_orig)
print [len(i) for i in ener_orig]
print len(ener_orig)
print [len(i) for i in dens_orig]
print len(dens_orig)
print [len(i) for i in H]
print len(H)
"""
"""

state_coords = params
param_types = ['epsilon','rmin_half']
all_xyz = xyz_orig
filename = 'packmol_boxes/cyclohexane_250.pdb'
pdb = PDBFile(filename)
# Re-evaluate energies from original parameter state to subsample and use that indexing for the other energies and MBAR calculations
for ii,value in enumerate(all_xyz): 
    print "Number of energy evaluations: %s" %(len(state_coords))
    print "starting energy evaluations"
    D = OrderedDict()
    for i,val in enumerate(state_coords):
        D['State' + ' ' + str(i)] = [["[#6X4:1]",param_types[j],val[j]] for j in range(len(param_types))]#len(state_orig))]
    D_mol = {'cyclohexane' : D}
  
    # Return energies from native parameter state to use as timeseries to subsample
    energies, u = new_param_energy(all_xyz[ii], D_mol, pdb.topology, vecs_orig[ii], T = 293.15)

    ts_sub,N_k_sub,xyz_sub,ind_subs = subsampletimeseries([energies],[all_xyz[ii]],[len(all_xyz[ii])])
        
ind_subs = ind_subs[0]
#vecs_sub = []
#for i,j in zip(vecs_orig,ind_subs):
#    print i
#    print j
#    vec_sub.append(i[j])
vecs_sub = [vecs_orig[0][j] for j in ind_subs]
vol_sub = [vol_orig[0][j] for j in ind_subs]

######################################################################################### WORKS^^
# Define new parameter states we wish to evaluate energies at
more = [['0.05','1.9080'],['0.1094','1.5'],['0.15','1.9080']]
for i in more:
    state_coords.append(i)


A_expectations = [[] for a in state_coords]
A_var_expectations = [[] for a in state_coords]
A_var_alt_expectations = [[] for a in state_coords]
dA_expectations = [[] for a in state_coords]
dA_var_expectations = [[] for a in state_coords]

for ii,value in enumerate(xyz_sub):
    MBAR_moves = state_coords
    print "Number of MBAR calculations: %s" %(len(MBAR_moves))
    print "starting MBAR calculations"
    D = OrderedDict()
    for i,val in enumerate(MBAR_moves):
        D['State' + ' ' + str(i)] = [["[#6X4:1]",param_types[j],val[j]] for j in range(len(param_types))]#len(state_orig))]
    D_mol = {'cyclohexane' : D} 
    print D_mol
    # Produce the u_kn matrix for MBAR based on the subsampled configurations
    E_kn, u_kn = new_param_energy(xyz_sub[ii], D_mol, pdb.topology, vecs_sub, T = 293.15)
    
    print u_kn
     
    K,N = np.shape(u_kn)
    print K,N
    N_k = np.zeros(K)
    N_k[0] = N

    # Initialize MBAR with Newton-Raphson
    # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
    initial_f_k = None # start from zero
    mbar = mb.MBAR(u_kn, N_k, verbose=False, relative_tolerance=1e-12,initial_f_k=initial_f_k)
    
    #for ind_E,E_k in enumerate(E_kn[ii]):
        #A_kn = [A_kn[ind_sub] for ind_sub in ind_subs[0]]
        #A_kn_var1 = [(A - np.average(A))**2 for A in A_kn]
        #------------------------------------------------------------------------
        # Compute Expectations for mass density and heat capacity (C_p)
        # Mass density (D) is easy, but recall that we may calculate C_p by thermal fluctuations:
        #
        #    C_p = (<E^2> - <E>^2)/(k_B * T^2)
        #
        #    or
        #
        #    C_p = <(E - <E>)^2>/(_B * T^2)
        #
        # We can add a number of harmonic quantum corrections, but we'll get to those later
        #------------------------------------------------------------------------
    (E_expect, dE_expect) = mbar.computeExpectations(E_kn,state_dependent = True)
    (E2_expect, dE2_expect) = mbar.computeExpectations(E_kn**2,state_dependent = True)
    (Vol_expect,dVol_expect) = mbar.computeExpectations(vol_sub,state_dependent = False)
    #(A_expect, dA_expect) = mbar.computeExpectations(A_kn,state_dependent = False)

    #A_kn_var = [(A - A_expect[0])**2 for A in A_kn]

    #(A_var_expect, dA_var_expect) = mbar.computeExpectations(A_kn_var,state_dependent = False)
    #(A_var_expect_alt,dA_var_expect_alt) = mbar.computeExpectations(A_kn_var1,state_dependent = False)

    #A_expectations[ii].append(A_expect)
    #A_var_expectations[ii].append(A_var_expect)
    #A_var_alt_expectations[ii].append(A_var_expect_alt)
    #A_expectations[ii].append(dA_expect)
    #dA_var_expectations[ii].append(dA_var_expect)
    #set_trace()
    E_fluc_input = np.array([(E_kn[i][:] - E_expect[i])**2 for i in range(len(E_expect))])
    (E_fluc_expect, dE_fluc_expect) = mbar.computeExpectations(E_fluc_input,state_dependent = True)
    T = 293.15 #K
    C_p_expect_meth1 = (E2_expect - E_expect**2)/(kB * T**2)
    C_p_expect_meth2 = E_fluc_expect/(kB * T**2)
"""  

files = nc.glob('*_cont.nc')
file_strings = [i.rsplit('_',1)[0] for i in files]
file_str = set(file_strings)

file_tups_traj = [[i+'.nc',i+'_cont.nc'] for i in file_str]

file_tups_sd = [[i+'.csv',i+'_cont.csv'] for i in file_str]
params = [i.rsplit('.',1)[0].rsplit('_') for i in files]
params = [[i[3][7:],i[5][4:]] for i in params]
MMcyc = 0.08416 #kg/mol

states_traj = [[] for i in file_tups_traj]
states_sd = [[] for i in file_tups_sd]
xyz_orig = [[] for i in file_tups_traj]
vol_orig = [[] for i in file_tups_traj]
lens_orig = [[] for i in file_tups_traj]
angs_orig = [[] for i in file_tups_traj]
ener_orig = [[] for i in file_tups_sd]
dens_orig = [[] for i in file_tups_sd]
vecs_orig = [[] for i in file_tups_sd]
burnin = 1800

print 'Analyzing Cyclohexane neat liquid trajectories'
for j,i in enumerate(file_tups_traj):
    #print 'Analyzing trajectory %s of %s'%(j+1,len(file_tups_traj)+1)
    for ii in i:
        try:
            data, xyz, lens, angs = read_traj(ii,burnin)
            #set_trace()
        except IndexError:
            print "The trajectory had fewer than %s frames" %(burnin)
            continue
        for m,n in zip(lens,angs):
            vecs = md.utils.lengths_and_angles_to_box_vectors(float(m[0]),float(m[1]),float(m[2]),float(n[0]),float(n[1]),float(n[2]))
            vecs_orig[j].append(vecs)
        for pos in xyz:
            xyz_orig[j].append(pos)
        for l in lens:
            vol_orig[j].append(np.prod(l)) #A**3
            lens_orig[j].append(l)
        for ang in angs:
            angs_orig[j].append(ang)
    states_traj[j].append(i[0].rsplit('.',1)[0])

print [np.average(i) for i in vol_orig]
#print Vol_expect
ener2_orig = [[] for i in file_tups_sd]
for j,i in enumerate(file_tups_sd):
    try:
        datasets = [pd.read_csv(ii,sep=',')[burnin:] for ii in i]
        merged = pd.concat(datasets)
    except IndexError:
        print "The state data record had fewer than 250 frames"
    for e in merged["Potential Energy (kJ/mole)"]:
        ener_orig[j].append(e)
        ener2_orig[j].append(e**2)
    for dens in merged["Density (g/mL)"]:
        dens_orig[j].append(dens)
    states_sd[j].append(i[0].rsplit('.',1)[0])

#print ener_orig

print [(np.average(j) - np.average(i)**2)/(kB * T**2) for j,i in zip(ener2_orig,ener_orig)]
print [np.average((i - np.average(i))**2)/(kB * T**2) for i in ener_orig]
print C_p_expect_meth1
print C_p_expect_meth2
