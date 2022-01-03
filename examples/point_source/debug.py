import matplotlib.pyplot as plt
import stompy.grid.unstructured_grid as ugrid
import xarray as xr
##

ds=xr.open_dataset('rundata/Estuary_SUNTANS.nc.nc.0')

##

g=ugrid.UnstructuredGrid.from_ugrid(ds)

##

plt.figure(1).clf()
fig,axs=plt.subplots(2,1,sharex=True,sharey=True,num=1)

# zoom=(432.3059924308611, 539.4068707185047, 399.26002512857804, 442.8560563918511)
zoom=(96.30677537145674, 485.4896920988444, 160.29477373056574, 314.58686985096375)

sel=dict(time=2100,Nk=8)

salt=ds.salt.isel(**sel).values
ccoll=g.plot_cells(values=salt,
                   cmap='jet',clim=[0.95,1.05],ax=axs[0])
                   #labeler=lambda i,r: "%d: %.3f"%(i,salt[i]),
#                   clip=zoom)

if 0:
    dzz=ds.dzz.isel(**sel).values
    ccoll2=g.plot_cells(values=dzz,
                        cmap='jet',clim=[0,0.001],ax=axs[1])
    plt.colorbar(ccoll2,ax=axs[1],label='dzz')
if 1:
    w=ds.w.isel(time=sel['time'],Nkw=sel['Nk']).values
    ccoll2=g.plot_cells(values=w,
                        cmap='PuOr',clim=[-0.01,0.01],ax=axs[1])
    plt.colorbar(ccoll2,ax=axs[1],label='w')



plt.colorbar(ccoll,ax=axs[0],label='salt')

axs[0].axis('equal')
axs[0].axis(zoom)

