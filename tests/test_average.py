"""
Tests related to average output
"""

import os
from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.suntans import sun_driver

import xarray as xr

try:
    sun_driver.SuntansModel.sun_bin_dir=os.path.join(os.path.dirname(__file__),"../main")
except NameError:
    sun_driver.SuntansModel.sun_bin_dir="../main"

# Not very clean... but make sure we're using system mpi.
sun_driver.SuntansModel.mpi_bin_dir="/usr/bin/"

def base_model_setup(L=100,dx=10):
    g=unstructured_grid.UnstructuredGrid(max_sides=4)
    g.add_rectilinear([0,0],[L,10],1+int(L/dx),2)

    g.add_cell_field('depth',-6*np.ones(g.Ncells()))

    model=sun_driver.SuntansModel()
    model.num_procs=1
    model.load_template('point_source_test.dat')
    model.set_grid(g)
    model.run_start=np.datetime64("2018-01-01 00:00")
    model.run_stop =np.datetime64("2018-01-02 00:00")

    # 50m2, 5.0 m3/s, 0.1m/s
    source=sun_driver.FlowBC(name='inflow',geom=np.array([ [0,0], [0,10]]),
                             Q=10.0)
    outlet=sun_driver.StageBC(name='outflow',geom=np.array([ [100,0], [100,10]]),
                              z=-1)
    model.add_bcs([source,outlet])
    
    model.set_run_dir('rundata_2d_average', mode='pristine')
    model.projection='EPSG:26910'

    model.config['dt']=2.5
    model.config['ntout']=100
    model.config['ntaverage']=100
    model.config['calcaverage']=1
    model.config['averageNetcdfFile']="average.nc"
    model.config['thetaramptime']=100
    
    model.config['Cmax']=2
    model.config['Nkmax']=1
    model.config['stairstep']=0
    model.config['mergeArrays']=0
    return model
    
def test_2d_average():
    """
    Create a channel with steady flow, compare to average output
    """
    model=base_model_setup()
    g=model.grid
    model.write()

    model.ic_ds.eta.values[:]=-1
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()

    # The actual tests
    avg=xr.open_dataset(model.avg_outputs()[0])
    ds=xr.open_dataset(model.map_outputs()[0])

    # There had been issues with strides of the grid->edges array.
    # This is just a sanity check on the grid data in the average file
    assert np.all(avg.edges.values>=0)
    # while it's unclear whether the average output should have the same
    # number or one fewer timesteps, this test takes a stand.
    assert np.all( avg.time.values == ds.time.values )

    # Calculate flow in 3 ways, make sure that at steady state they
    # agree.
    
    # avoid the very end since it will have zero velocity (stage BC).
    j=g.select_edges_nearest([40,5]) # 11
    c=g.edge_to_cells(j)

    # Average flow
    Q_U_F=(avg.n1*avg.U_F).isel(Ne=j,Nk=0).values
    # Average depth, width and velocity
    w=avg.df.isel(Ne=j).values
    Q_avg_uc=w*((avg.dv+avg.eta)*avg.uc).isel(Nk=0,Nc=c[0]).values
    # Instantaneous depth and velocity
    Q_uc=w*((ds.dv+ds.eta)*ds.uc).isel(Nk=0,Nc=c[0]).values 

    # with the faster ramp up time, these are not quite the same, so
    # grab end of the time series where it's steady.
    assert np.abs(Q_U_F - Q_avg_uc)[-20].max() < 1e-4
    # compare last time step when nearly steady
    assert (Q_uc - Q_avg_uc)[-1] / Q_uc[-1] < 1e-4

    if 0:
        import matplotlib.pyplot as plt
        plt.figure(3).clf()
        plt.plot(avg.time,Q_U_F,label="U_F",lw=3)
        plt.plot(avg.time,Q_avg_uc,lw=5, alpha=0.5, label="uc")
        plt.plot(ds.time,Q_uc, label="map-uc")
        plt.legend()

    return model


def test_2d_average_merged():
    """
    Similar setup to test_2d_average, but check topology 
    in the merged grid
    """
    model=base_model_setup()
    g=model.grid
    model.num_procs=4
    model.run_stop =np.datetime64("2018-01-01 02:00")
    
    model.config['mergeArrays']=1
    model.write()

    model.ic_ds.eta.values[:]=-1
    model.write_ic_ds()

    model.partition()
    model.sun_verbose_flag='-v'
    model.run_simulation()

    # The actual tests
    avg=xr.open_dataset(model.avg_outputs()[0])
    ds=xr.open_dataset(model.map_outputs()[0])

    # There had been issues with strides of the grid->edges array.
    # This is just a sanity check on the grid data in the average file
    assert np.all(avg.edges.values>=0)
    # while it's unclear whether the average output should have the same
    # number or one fewer timesteps, this test takes a stand.
    assert np.all( avg.time.values == ds.time.values )

    # Calculate flow in 3 ways, make sure that at steady state they
    # agree.
    
    # avoid the very end since it will have zero velocity (stage BC).
    j=g.select_edges_nearest([40,5]) # 11
    c=g.edge_to_cells(j)

    # Average flow
    Q_U_F=(avg.n1*avg.U_F).isel(Ne=j,Nk=0).values
    # Average depth, width and velocity
    w=avg.df.isel(Ne=j).values
    Q_avg_uc=w*((avg.dv+avg.eta)*avg.uc).isel(Nk=0,Nc=c[0]).values
    # Instantaneous depth and velocity
    Q_uc=w*((ds.dv+ds.eta)*ds.uc).isel(Nk=0,Nc=c[0]).values 

    # with the faster ramp up time, these are not quite the same, so
    # grab end of the time series where it's steady.
    assert np.abs(Q_U_F - Q_avg_uc)[-20].max() < 1e-4
    # compare last time step when nearly steady
    assert (Q_uc - Q_avg_uc)[-1] / Q_uc[-1] < 1e-4

    if 0:
        import matplotlib.pyplot as plt
        plt.figure(3).clf()
        plt.plot(avg.time,Q_U_F,label="U_F",lw=3)
        plt.plot(avg.time,Q_avg_uc,lw=5, alpha=0.5, label="uc")
        plt.plot(ds.time,Q_uc, label="map-uc")
        plt.legend()

    # might as well test ability to load model back in.
    model=sun_driver.SuntansModel.load('rundata_2d_average')

    avg_ds=xr.open_dataset(model.avg_outputs()[0])
    map_ds=xr.open_dataset(model.map_outputs()[0])
    for ds in [avg_ds,map_ds]:
        ds.load() 
        ds.close() # avoids some pitfalls with repeat runs

    avg_g=unstructured_grid.UnstructuredGrid.from_ugrid(avg_ds)
    map_g=unstructured_grid.UnstructuredGrid.from_ugrid(map_ds)

    # some topology checks on the merged grids
    for ds,g in [(avg_ds,avg_g),
                 (map_ds,map_g)]:
        assert np.allclose(ds['xp'].values,g.nodes['x'][:,0])
        assert np.allclose(ds['yp'].values,g.nodes['x'][:,1])
        assert np.all( ds['cells']==g.cells['nodes'] )

        c2e=np.array( [g.cell_to_edges(c) for c in range(g.Ncells())] )
        assert np.all( ds['face'].values==c2e)
        assert np.all( ds['edges'].values==g.edges['nodes'] )

        # Frustrating, but unstructured_grid and suntans have opposite sign
        # conventions for edge normals.  at least make sure it's consistent
        assert np.all( ds['n1'].values == -g.edges_normals()[:,0] )
        assert np.all( ds['n2'].values == -g.edges_normals()[:,1] )

        # Those all pass
        e2c=g.edge_to_cells()
        e2c[e2c<0]=-1 # standardize the missing values
        assert np.all(ds.grad.values==e2c)



