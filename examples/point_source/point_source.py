"""
Development and testing of point sources, i.e. discharge
located at bed
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

n1=g.select_nodes_nearest([300,0])
n2=g.select_nodes_nearest([300,500])

model=sun_driver.SuntansModel()
model.load_template('sun-template.dat')
model.set_grid(g)

source=dfm.SourceSinkBC(name='inflow',geom=np.array([200,200]),
                        z=-10,Q=1.0)
model.add_bcs(source)

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

#

# manually write the point source
Npoint=0
point_cell=[]
point_layer=[]
point_Q=[]
point_S=[]
point_T=[]
Nt=model.bc_ds.dims['Nt']

for key in model.bc_point_sources.keys():
    (c,k)=key
    print("Point source for cell=%d, k=%d"%(c,k))
    assert 'Q' in model.bc_point_sources[key]

    # this will be handled by combine_items
    das=model.bc_point_sources[key]['Q']
    assert len(das)==1
    da=das[0]
    
    # when moving this into the core code, can re-use the
    # time-handling utilities there.
    # this is interp_time(...)
    assert 'time' not in da.dims,"wait for that"
    data=da.values * np.ones( (Nt,)+da.values.shape )
    point_cell.append(c)
    point_layer.append(k)
    point_Q.append(data) 
    print("punting on point source salinity and temp")
    point_S.append(0*np.ones_like(data))
    point_T.append(4*np.ones_like(data))

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

