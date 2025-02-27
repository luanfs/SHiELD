import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
from matplotlib.colors import ListedColormap



#----------------------------------------------------------------------------------------------
#source https://unidata.github.io/python-gallery/examples/Precipitation_Map.html
cmap_data = [(1.0, 1.0, 1.0),
             (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.250980406999588, 0.250980406999588, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0)]
 
# Create the colormap using ListedColormap
cmap_precp = ListedColormap(cmap_data)
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield(q, title, filename, filepath, colormap, qmin, qmax):
    print('plotting '+filename)
    # Figure quality
    dpi = 100

    # Figure format
    figformat='png'

    # Map projection
    map_projection = "mercator"
    #map_projection = "sphere"

    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1600/dpi, 1000/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()
    elif map_projection == "sphere":
        plt.figure(figsize=(1000/dpi,1000/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)

    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 19, 'color': 'black'}
    gl.ylabel_style = {'size': 19, 'color': 'black'}

 
    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    N = np.shape(q)[0]
    # plot for each tile
    for tile in range(0,6):
        # Grid file to be opened
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the file
        grid = xr.open_dataset(grid_file , decode_times=False)

        # Get grid
        lon = grid['grid_lon']
        lat = grid['grid_lat']

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.2
        plt.rcParams["axes.axisbelow"] = True

        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-', linewidth=lw, color=colors[tile], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,tile], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)

    plt.title(title,  fontsize=19)
    #print(qmin,qmax)
    # Plot colorbar
    cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.046, pad=0.04, format='%.1f')
    #cax,kw = colorbar.make_axes(ax,orientation='horizontal' , fraction=0.09, pad=0.04, shrink=0.9, format='%.0e')

    cb=plt.colorbar(im, cax=cax, extend='both',**kw)

    ticks = np.linspace(qmin, qmax, num=5)
    cb.set_ticks(ticks)
    cb.ax.tick_params(labelsize=22)
    plt.savefig(filename+'.'+figformat, format=figformat)
    plt.show()
    plt.close()
    #exit()
#-----------------------------------------------------------------------------------------
