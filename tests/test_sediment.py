"""
Tests related to point sources
"""

import os
from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.suntans import sun_driver

import xarray as xr

sun_driver.SuntansModel.sun_bin_dir=os.path.join(os.path.dirname(__file__),"../main")

def base_model():
    """
    A channel, uniform bathy.
    """
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[1000,100],21,3)
    g.add_cell_field('depth',-6*np.ones(g.Ncells()))

    model=sun_driver.SuntansModel()
    model.load_template('point_source_test.dat')
    model.set_grid(g)
    model.run_start=np.datetime64("2018-01-01 00:00")
    model.run_stop =np.datetime64("2018-01-01 10:00")

    dt=np.timedelta64(600,'s')
    times=np.arange(model.run_start-dt,model.run_stop+2*dt,dt)
    secs=(times-times[0])/np.timedelta64(1,'s')
    eta0=-2
    eta_bc=sun_driver.StageBC(name="eta_bc",
                              geom=np.array([ [0,0],
                                              [0,100]]),
                              z=eta0)

    model.add_bcs(eta_bc)
    model.add_bcs( [sun_driver.ScalarBC(parent=eta_bc,scalar="S",value=1),
                    sun_driver.ScalarBC(parent=eta_bc,scalar="T",value=1)] )

    model.set_run_dir('rundata_sedi_3d', mode='pristine')
    model.projection='EPSG:26910'
    model.num_procs=1 # test single first
    model.config['dt']=30
    model.config['ntout']=5
    model.config['Cmax']=30
    model.config['Nkmax']=10
    model.config['stairstep']=0
    model.config['mergeArrays']=0

    # Sediment:
    model.config['computeSediments']=1   # whether include sediment model
    
    model.config['Nlayer']=0 # Number of bed layers (MAX = 5)
    model.config['Nsize']=3  # Number of sediment fractions (Max = 3)
    model.config['TBMAX']=1  # whether to output tb for each cell
    model.config['SETsediment']=0       # When Nlayer>5 or Nsize>3, SETsediment=1 to use SetSediment 
    model.config['WSconstant']=1       # if 1, ws for sediment = ws0
    model.config['readSediment']=0       # if 1, read sediment concentration data as IC. only work with Nsize==1
    model.config['bedInterval']=150     # the interval steps to update bed change
    model.config['bedComplex']=0       # whether to consider the possibility to flush away a whole layer
    model.config['ParabolKappa']=0       # whether to use parabolic tubulent diffusivity
    model.config['Ds90']=0.000008    # ds90 for calculation of erosion taub
    model.config['Ds1']=0.00000057   # sediment diameter for fraction No.1 (m)
    model.config['Ds2']=0.0002    # sediment diameter for fraction No.2                   
    model.config['Ds3']=0.0002    # sediment diameter for fraction No.3
    model.config['Ws01']= 0.0001     # constant settling velocity for fraction No.1 (m/s)
    model.config['Ws02']= 0.01     # constant settling velocity for fraction No.2
    model.config['Ws03']=-0.01     # constant settling velocity for fraction No.3
    model.config['Gsedi1']=2.65      # relative density for fraction No.1
    model.config['Gsedi2']=2.65      # relative density for fraction No.2
    model.config['Gsedi3']=2.65      # relative density for fraction No.3
    model.config['Prt1']=1       # Prandtl Number for fraction No.1
    model.config['Prt2']=1       # Prandtl Number for fraction No.2
    model.config['Prt3']=1       # Prandtl Number for fraction No.3
    model.config['Consolid1']=0.00    # Consolidation rate (g/m^2/s) for layer No.1
    model.config['E01']=0.1      # Basic Erosion Rate Constant (g/m^2/s) for layer No.1
    model.config['Taue1']=0.0       # Erosion Critical Shear Stress (N/m^2) for layer No.1
    # Taud1 isn't actually used in current code.  Deposition always occurs at the rate
    # of w_s.
    model.config['Taud1']=0.0       # Deposition Critical Shear Stress (N/m^2) for layer No.1
    model.config['Drydensity1']=530000    # Dry density (g/m^3) for layer No.1
    model.config['Thickness1']=0.0      # Thickness (m) for layer No.1
    model.config['softhard1']=1         # 0 soft or hard for layer No.1 to decide how to calculate erosion
    model.config['Bedmudratio1']=0.8       # Bed mud ratio for layer No.1
    model.config['Chind']=1000000   # Concentration (in volumetric fraction) criterion for hindered settling velocity
    model.config['Cfloc']=500000    # Concentration (in volumetric fraction) criterion for flocculated settling velcoity
    model.config['k']=0.0002    # Constant coefficient for settling velocity as a function of conc.

    model.config['Sediment1File']='sedi1.dat'
    model.config['Sediment2File']='sedi2.dat'
    model.config['Sediment3File']='sedi3.dat'
    model.config['LayerFile']='sedilayer.dat'
    model.config['tbFile']='Seditb.dat'
    model.config['tbmaxFile']='Seditbmax.dat'
    return model

def test_sed_quiescent():
    """
    quiescent, uniform ic, sed1 settles down, sed2 settles up.
    """
    model=base_model()
    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=model.bcs[0].z
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0

    # let sed1 default to 0
    for sed in ['sed2','sed3']:
        model.ic_ds[sed]=model.ic_ds['salt'].copy()
        model.ic_ds[sed].values[:]=100.0
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model

def test_sed_ic():
    """
    spatially variable IC.
    """
    model=base_model()
    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=model.bcs[0].z
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0

    # let sed1 default to 0
    for sed in ['sed2','sed3']:
        model.ic_ds[sed]=model.ic_ds['salt'].copy()
        if sed=='sed2':
            sed_grad=model.ic_ds.xv.values / 1000.
        else:
            sed_grad=(1000-model.ic_ds.xv.values) / 1000.
        # time,Nk,Nc
        model.ic_ds[sed].values[:,:,:]=sed_grad[None,None,:]
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model

def test_sed_bc():
    """
    basic BC
    """
    model=base_model()

    # slow down the settling velocities to get some advection mixed in
    model.config['Ws01']= 0.0001     # constant settling velocity for fraction No.1 (m/s)
    model.config['Ws02']= 0.0001     # constant settling velocity for fraction No.2
    model.config['Ws03']=-0.0001     # constant settling velocity for fraction No.3

    geom=np.array([ [1000,0],
                    [1000,100]])
    Q_bc=sun_driver.FlowBC(name="Q_bc",
                           geom=geom,
                           Q=100.0)
    
    model.add_bcs( [Q_bc,
                    sun_driver.ScalarBC(parent=Q_bc,scalar="S",value=2),
                    sun_driver.ScalarBC(parent=Q_bc,scalar="T",value=0)] )

    point_bc=sun_driver.SourceSinkBC(name="bedPoint",
                                     geom=np.array([500,20]),
                                     Q=10)
    model.add_bcs( [point_bc,
                    sun_driver.ScalarBC(parent=point_bc,scalar="S",value=3),
                    sun_driver.ScalarBC(parent=point_bc,scalar="T",value=3)] )
                    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=model.bcs[0].z
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0

    model.write_ic_ds()
        
    for sed_idx in range(3):
        name="sed%d"%(sed_idx+1)
        model.bc_ds['boundary_'+name]=model.bc_ds['boundary_S'].copy()
        model.bc_ds[name]=model.bc_ds['S'].copy()
        model.bc_ds['point_'+name]=model.bc_ds['point_S'].copy()
        
        model.bc_ds['boundary_'+name].values[:]=sed_idx*1.0
        model.bc_ds[name].values[:]=0.0
        model.bc_ds['point_'+name].values[:]=sed_idx*2.0

    model.write_bc_ds()
    
    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model

def test_sed_mpi():
    model=base_model()
    
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[1000,100],101,11)
    g.add_cell_field('depth',-6*np.ones(g.Ncells()))

    model.set_grid(g)

    model.num_procs=4
    model.config['dt']=3.

    # slow down the settling velocities to get some advection mixed in
    model.config['Ws01']= 0.0001     # constant settling velocity for fraction No.1 (m/s)
    model.config['Ws02']= 0.0005     # constant settling velocity for fraction No.2
    model.config['Ws03']=-0.0005     # constant settling velocity for fraction No.3

    geom=np.array([ [1000,0],
                    [1000,100]])
    Q_bc=sun_driver.FlowBC(name="Q_bc",
                           geom=geom,
                           Q=100.0)
    
    model.add_bcs( [Q_bc,
                    sun_driver.ScalarBC(parent=Q_bc,scalar="S",value=2),
                    sun_driver.ScalarBC(parent=Q_bc,scalar="T",value=0)] )

    point_bc=sun_driver.SourceSinkBC(name="bedPoint",
                                     geom=np.array([500,20]),
                                     Q=10)
    model.add_bcs( [point_bc,
                    sun_driver.ScalarBC(parent=point_bc,scalar="S",value=3),
                    sun_driver.ScalarBC(parent=point_bc,scalar="T",value=3)] )
                    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=model.bcs[0].z
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0

    model.write_ic_ds()
        
    for sed_idx in range(3):
        name="sed%d"%(sed_idx+1)
        model.bc_ds['boundary_'+name]=model.bc_ds['boundary_S'].copy()
        model.bc_ds[name]=model.bc_ds['S'].copy()
        model.bc_ds['point_'+name]=model.bc_ds['point_S'].copy()
        
        model.bc_ds['boundary_'+name].values[:]=sed_idx*1.0
        model.bc_ds[name].values[:]=0.0
        model.bc_ds['point_'+name].values[:]=sed_idx*2.0

    model.write_bc_ds()
    
    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model
    

model=test_sed_mpi()
