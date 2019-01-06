"""
Development and testing of independently-specified
edge depths
"""

from stompy.grid import unstructured_grid
import numpy as np
from stompy.model.delft import dflow_model
from stompy.model.suntans import sun_driver

import six
six.moves.reload_module(dflow_model)
six.moves.reload_module(sun_driver)


g=unstructured_grid.UnstructuredGrid(max_sides=4)

ret=g.add_rectilinear([0,0],[500,500],50,50)

g.add_cell_field('depth',-10*np.ones(g.Ncells()))
g.add_edge_field('edge_depth',-10*np.ones(g.Nedges()))


n1=g.select_nodes_nearest([300,0])
n2=g.select_nodes_nearest([300,500])

edges=g.shortest_path( n1,n2, return_type='edges' )
# sloping weir, initially totally dry, will overtop
# during run
# test both signs of the edge depth
g.edges['edge_depth'][edges]= np.linspace(-2,2,len(edges))

model=sun_driver.SuntansModel()
model.use_edge_depths=True
model.load_template('sun-template.dat')
model.set_grid(g)

model.num_procs=4

inflow=sun_driver.FlowBC(name='inflow',
                         geom=np.array([ [0,0],[0,100]]),
                         Q=50.0)
model.add_bcs(inflow)

model.set_run_dir('rundata', mode='pristine')
model.run_start=np.datetime64("2018-01-01 00:00")
model.run_stop =np.datetime64("2018-01-01 20:00")
model.projection='EPSG:26910'
model.sun_bin_dir="/home/rusty/src/suntans/main"

# for 2D, dt=30 is okay.
# for 3D, start getting some CmaxW, so scale it back.
#     10 or 15 would probably be okay. 30 is unstable.
model.config['dt']=5
model.config['Cmax']=30
model.config['Nkmax']=10
model.config['stairstep']=0

model.write()

model.ic_ds.eta.values[:]=-2.1
model.write_ic_ds()

model.partition()
model.run_simulation()

