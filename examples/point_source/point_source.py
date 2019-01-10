"""
Development and testing of point sources, i.e. discharge
located at bed.
"""

from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.delft import dflow_model
from stompy.model.suntans import sun_driver

import six
six.moves.reload_module(dflow_model)
six.moves.reload_module(sun_driver)

# Create square grid
g=unstructured_grid.UnstructuredGrid(max_sides=4)
g.add_rectilinear([0,0],[500,500],50,50)

g.add_cell_field('depth',-10*np.ones(g.Ncells()))

n1=g.select_nodes_nearest([300,0])
n2=g.select_nodes_nearest([300,500])

model=sun_driver.SuntansModel()
model.load_template('sun-template.dat')
model.set_grid(g)

source=dfm.SourceSinkBC(name='inflow',geom=np.array([200,200]),
                        z=-10,Q=1.0)

model.add_bcs(source)
# HERE - this doesn't feel like a reasonable way to handle
# scalars at source/sinks. it's not enough, though, to just specify
# concentration at a location, because there could be multiple sources
# mapping to the same cell/layer.
# in that sense, it would be better to do something like this, but then
# change flow bcs to do the same.  I guess it is the difference between
# adding and setting?  It's going to get ugly -- back off and just
# set the concentrations, but at least allow reusing the source info.
model.add_bcs( [source.scalar_bc(salt=0),
                source.scalar_bc(temperature=4)] )

model.set_run_dir('rundata', mode='pristine')
model.run_start=np.datetime64("2018-01-01 00:00")

model.run_stop =np.datetime64("2018-01-01 01:00")
model.projection='EPSG:26910'
model.sun_bin_dir="/home/rusty/src/suntans/main"

model.config['dt']=10
model.config['Cmax']=30
model.config['Nkmax']=5
model.config['stairstep']=0

model.write()

# easier to debug when the top layer is wet
model.ic_ds.eta.values[:]=-0.5
model.ic_ds.salt.values[:]=1.0
model.ic_ds.temp.values[:]=1.0
model.write_ic_ds()

bc_ds=model.bc_ds
bc_ds['point_cell']=('Npoint',), point_cell
bc_ds['point_layer']=('Npoint',), point_layer
bc_ds['point_Q']=('Nt','Npoint'), np.stack(point_Q,axis=-1)
bc_ds['point_S']=('Nt','Npoint'), np.stack(point_S,axis=-1)
bc_ds['point_T']=('Nt','Npoint'), np.stack(point_T,axis=-1)

model.write_bc_ds()

#
model.partition()
model.sun_verbose_flag='-v'
model.run_simulation()

