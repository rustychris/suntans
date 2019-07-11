"""
Tests related to nudging of temperature, evaporation, rain.
"""

import os
from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.suntans import sun_driver

import xarray as xr

sun_driver.SuntansModel.sun_bin_dir=os.path.join(os.path.dirname(__file__),"../main")

def base_model():
    """
    Channel
    """
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[1000,100],21,3)
    
    g.add_cell_field('depth',-6*np.ones(g.Ncells()))

    model=sun_driver.SuntansModel()
    model.load_template('point_source_test.dat')
    model.set_grid(g)
    model.run_start=np.datetime64("2018-01-01 00:00")
    model.run_stop =np.datetime64("2018-01-05 00:00")

    dt=np.timedelta64(600,'s')
    times=np.arange(model.run_start-dt,model.run_stop+2*dt,dt)
    secs=(times-times[0])/np.timedelta64(1,'s')
    eta0=-1
    eta_bc=sun_driver.StageBC(name="eta_bc",
                              geom=np.array([ [0,0],
                                              [0,100]]),
                              z=eta0)

    model.add_bcs(eta_bc)
    model.add_bcs( [sun_driver.ScalarBC(parent=eta_bc,scalar="S",value=34),
                    sun_driver.ScalarBC(parent=eta_bc,scalar="T",value=15)] )

    model.set_run_dir('rundata_met_3d', mode='pristine')
    model.projection='EPSG:26910'
    model.num_procs=1 # test single first
    model.config['dt']=120
    model.config['ntout']=5
    model.config['Cmax']=30
    model.config['Nkmax']=10
    model.config['stairstep']=0
    model.config['mergeArrays']=0

    return model

def test_met_quiescent():
    """
    quiescent, uniform ic
    """
    model=base_model()

    model.config['metmodel']=1

    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=model.bcs[0].z
    model.ic_ds.salt.values[:]=34.0
    model.ic_ds.temp.values[:]=15.0
    
    model.write_ic_ds()

    model.met_ds['Tair'].values[:]=20
    model.met_ds['cloud'].values[:]=0.5
    model.met_ds['RH'].values[:]=100.0
    model.write_met_ds()
    
    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model


def test_met_nudge():
    """
    test metmodel=5, nudge to Tair
    This runs, has some minor w oscillations, but otherwise
    behaves.
    """
    model=base_model()
    model.run_start=np.datetime64("2018-06-10 00:00")
    model.run_stop =np.datetime64("2018-06-20 00:00")

    model.config['metmodel']=5
    model.config['rstretch']=1.1

    # quiescent first:
    model.bcs=[]
    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=-1.0
    model.ic_ds.salt.values[:]=34.0
    model.ic_ds.temp.values[:]=10.0
    
    model.write_ic_ds()

    model.met_ds['Tair'].values[:]=20
    model.write_met_ds()
    
    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()
    return model


def test_met_evap():
    """
    test metmodel=5, disable nudge to Tair, include evaporation
    """
    model=base_model()
    model.run_start=np.datetime64("2018-06-10 00:00")
    model.run_stop =np.datetime64("2018-06-20 00:00")

    model.config['metmodel']=5
    model.config['rstretch']=1.1
    model.config['Lsw']=0.0 # disable nudging

    # quiescent first:
    model.bcs=[]
    
    model.write()

    # leave some dry layers at the surface
    model.ic_ds.eta.values[:]=-1.0
    model.ic_ds.salt.values[:]=34.0
    model.ic_ds.temp.values[:]=10.0
    
    model.write_ic_ds()

    # Add data point for Tair and rain.
    model.met_ds['Tair'].values[:]=20
    # 5cm/day ~ 5e-7 m/s
    model.met_ds['rain'].values[:]=-1e-6
    
    model.write_met_ds()
    
    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()

    # did it actually work?
    ds=xr.open_dataset(model.map_outputs()[0])
    ds_last=ds.isel(time=-1)

    # make sure actual evaporation got reported, and
    # that the impose evaporation (negative rain) gets reported
    # here as positive EP, since EP == evaporation - precip.
    assert np.any( ds_last.EP.values > 0.0 )
    # salinity should be increase everywhere
    wet=ds_last.dzz.values>0
    assert np.all( ds_last.salt.values[ wet ] > 34.00 )
    assert np.abs( ds_last.temp.values[wet] - 10.0 ).max() < 0.01
    
    return model

model=test_met_evap()



