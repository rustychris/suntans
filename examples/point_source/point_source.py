"""
Development and testing of point sources, i.e. discharge
located at bed.

Also test 3D scalar to check for regressions.
"""

from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.delft import dflow_model
from stompy.model.suntans import sun_driver

import matplotlib.pyplot as plt
import xarray as xr

import six
six.moves.reload_module(dflow_model)
six.moves.reload_module(sun_driver)

# Create square grid
g=unstructured_grid.UnstructuredGrid(max_sides=4)
g.add_rectilinear([0,0],[500,500],50,50)

# wavy bed
half_wave=150
cc=g.cells_center()
cell_depth=-6 + np.cos(cc[:,0]*np.pi/half_wave) * np.cos(cc[:,1]*np.pi/half_wave)

g.add_cell_field('depth',cell_depth)

#plt.clf()
#g.plot_cells(values=cell_depth,cmap='jet')
#plt.axis('equal')

##


n1=g.select_nodes_nearest([300,0])
n2=g.select_nodes_nearest([300,500])

model=sun_driver.SuntansModel()
model.load_template('sun-template.dat')
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
# in that sense, it would be better to do something like this, but then
# change flow bcs to do the same.  I guess it is the difference between
# adding and setting?  It's going to get ugly -- back off and just
# set the concentrations, but at least allow reusing the source info.
model.add_bcs( [sun_driver.ScalarBC(parent=hill_source,scalar="S",value=1),
                sun_driver.ScalarBC(parent=hill_source,scalar="T",value=1)] )
model.add_bcs( [sun_driver.ScalarBC(parent=saddle_source,scalar="S",value=1),
                sun_driver.ScalarBC(parent=saddle_source,scalar="T",value=1)] )

model.set_run_dir('rundata', mode='pristine')
model.projection='EPSG:26910'
model.sun_bin_dir="/home/rusty/src/suntans/main"

model.config['dt']=2.5
model.config['ntout']=5
model.config['Cmax']=30
model.config['Nkmax']=10
model.config['stairstep']=0
model.config['mergeArrays']=0

model.write()

# easier to debug when the top layer is wet
model.ic_ds.eta.values[:]=eta_da.values[1]
model.ic_ds.salt.values[:]=1.0
model.ic_ds.temp.values[:]=1.0
model.write_ic_ds()

#
model.partition()
model.sun_verbose_flag='-v'
model.run_simulation()

