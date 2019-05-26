"""
Tests related to point sources
"""

import os
from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.suntans import sun_driver

import xarray as xr

sun_driver.SuntansModel.sun_bin_dir=os.path.join(os.path.dirname(__file__),"../main")

def test_point_source_1d():
    """
    Create a 1D water column, with a point source at the bed.
    Verifies that all cells with valid thickness (dzz>0.001m)
    have near unity salinity.
    """
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[10,10],2,2)

    g.add_cell_field('depth',-6*np.ones(g.Ncells()))

    model=sun_driver.SuntansModel()
    model.load_template('point_source_test.dat')
    model.set_grid(g)
    model.run_start=np.datetime64("2018-01-01 00:00")
    model.run_stop =np.datetime64("2018-01-01 10:00")

    source=sun_driver.SourceSinkBC(name='inflow',geom=np.array([5,5]),
                                   z=-10,Q=1.0)
    model.add_bcs(source)
    model.add_bcs( [sun_driver.ScalarBC(parent=source,scalar="S",value=1),
                    sun_driver.ScalarBC(parent=source,scalar="T",value=1)] )

    model.set_run_dir('rundata_point_1d', mode='pristine')
    model.projection='EPSG:26910'

    model.config['dt']=2.5
    model.config['ntout']=1
    model.config['Cmax']=30
    model.config['Nkmax']=10
    model.config['stairstep']=0
    model.config['mergeArrays']=0

    model.write()

    model.ic_ds.eta.values[:]=-5.999
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()

    ds=xr.open_dataset(model.map_outputs()[0])

    dzz=ds.dzz.values
    dzz_thresh=0.001
    salt=ds.salt.values
    salt_errors=salt[dzz>=dzz_thresh] - 1.0
    
    assert np.max( np.abs(salt_errors) )<0.001

    
def test_point_source_3d():
    """
    A square 3D grid with wetting and drying.
    the bathymetry is a wavy cos(x)*cos(y) function.
    two point sources and an oscillating freesurface.
    one point source at a saddle point, and the other
    on top of a bump which goes dry.

    This currently does pretty well, but it could be
    better -- see the test thresholds near the end which
    have been relaxed somewhat.
    """
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[500,500],50,50)

    # wavy bed
    half_wave=150
    cc=g.cells_center()
    cell_depth=-6 + np.cos(cc[:,0]*np.pi/half_wave) * np.cos(cc[:,1]*np.pi/half_wave)
    g.add_cell_field('depth',cell_depth)

    model=sun_driver.SuntansModel()
    model.load_template('point_source_test.dat')
    model.set_grid(g)
    model.run_start=np.datetime64("2018-01-01 00:00")
    model.run_stop =np.datetime64("2018-01-01 10:00")

    dt=np.timedelta64(600,'s')
    times=np.arange(model.run_start-dt,model.run_stop+2*dt,dt)
    secs=(times-times[0])/np.timedelta64(1,'s')
    eta_values=-6 + 0.75*np.cos(secs*np.pi/7200)
    eta_da=xr.DataArray(eta_values,coords=[ ('time',times) ])
    eta_bc=sun_driver.StageBC(name="eta_bc",geom=np.array([ [0,0],
                                                            [0,500]]),
                              z=eta_da)

    model.add_bcs(eta_bc)
    model.add_bcs( [sun_driver.ScalarBC(parent=eta_bc,scalar="S",value=1),
                    sun_driver.ScalarBC(parent=eta_bc,scalar="T",value=1)] )

    # point source that is typically dry
    hill_source=sun_driver.SourceSinkBC(name='inflow',geom=np.array([150,150]),
                                        z=-10,Q=1.0)

    # on a saddle
    saddle_source=sun_driver.SourceSinkBC(name='inflow',geom=np.array([220,225]),
                                          z=-10,Q=1.0)

    model.add_bcs(hill_source)
    model.add_bcs(saddle_source)
    model.add_bcs( [sun_driver.ScalarBC(parent=hill_source,scalar="S",value=1),
                    sun_driver.ScalarBC(parent=hill_source,scalar="T",value=1)] )
    model.add_bcs( [sun_driver.ScalarBC(parent=saddle_source,scalar="S",value=1),
                    sun_driver.ScalarBC(parent=saddle_source,scalar="T",value=1)] )

    model.set_run_dir('rundata_point_3d', mode='pristine')
    model.projection='EPSG:26910'
    model.num_procs=4
    model.config['dt']=2.5
    model.config['ntout']=5
    model.config['Cmax']=30
    model.config['Nkmax']=10
    model.config['stairstep']=0
    model.config['mergeArrays']=0

    model.write()

    model.ic_ds.eta.values[:]=eta_da.values[1]
    model.ic_ds.salt.values[:]=1.0
    model.ic_ds.temp.values[:]=1.0
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()

    for map_fn in model.map_outputs():
        ds=xr.open_dataset(map_fn)

        dzz=ds.dzz.values
        # This is not as strict as it should be.
        # should be 0.001 or less.
        dzz_thresh=0.05
        salt=ds.salt.values
        valid=np.isfinite(dzz) & (dzz>=dzz_thresh)
        salt_errors=salt[valid] - 1.0

        # This is not as strict as it should be.
        # we should be able to get this down to 1e-4 or
        # better.
        assert np.max( np.abs(salt_errors) )<0.01
        ds.close()

        
