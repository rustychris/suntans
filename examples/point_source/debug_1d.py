import matplotlib.pyplot as plt
import stompy.grid.unstructured_grid as ugrid
import xarray as xr
##

ds=xr.open_dataset('rundata_1d/Estuary_SUNTANS.nc.nc.0')

##

g=ugrid.UnstructuredGrid.from_ugrid(ds)

##

plt.figure(1).clf()
fig,axs=plt.subplots(3,1,sharex=True,sharey=True,num=1)

img1=axs[0].imshow( ds.salt.values[:,:,0].T, aspect='auto',cmap='jet',clim=[0.99,1.01],
                    interpolation='none')

plt.colorbar(img1,ax=axs[0])



img2=axs[1].imshow( ds.dzz.values[:,:,0].T, aspect='auto',cmap='jet',clim=[0,0.2],
                    interpolation='none')
plt.colorbar(img2,ax=axs[1])


img3=axs[2].imshow( ds.w.values[:,:,0].T, aspect='auto',cmap='seismic',clim=[-0.05,0.05],
                    interpolation='none')

plt.colorbar(img3,ax=axs[2])

axs[0].axis( (-13.710330434775543,
              305.6093689651833,
              10.142030945230957,
              -1.645655402704925) )
