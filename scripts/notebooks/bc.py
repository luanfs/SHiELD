import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
from matplotlib.colors import ListedColormap
import os
from plot_cs import *


#----------------------------------------------------------------------------------------------
# directories
hydro="nh"
inputdir='/gpfs/f5/scratch/Luan.Santos/gfdl_w/solo_'+hydro+'/'
outdir = '/gpfs/f5/scratch/Luan.Santos/gfdl_w/graphs_solo_'+hydro+'/'
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
#simulation parameters
N=96             # N
gtype = 0        # grid type
hord  = 8        # PPM option
modes = (32, 64) # single (32) and double (64) precisions"
Tf    = 15        # time of integration (days)
testname='bc'
basename= "C"+str(N)+"."+hydro+"."+testname+".g"+str(gtype)+".hord"+str(hord)+".T"+str(Tf)+"."
#----------------------------------------------------------------------------------------------

# Create filepaths list
filepaths = []
for mode in modes: 
    filepath = inputdir+basename+str(mode)+"bit"+"/rundir/"
    if( not os.path.exists(filepath) ):
       print("ERROR: filepath below does not exist")
       print(filepath)
       exit()
    filepaths.append(filepath)
    print(filepath)

# Get the number of plots
atmos_file = filepaths[0]+"atmos_daily.tile1.nc"
data = xr.open_dataset(atmos_file, decode_times=False)
times = data.time.values
nplots = len(times)

#-----------------------------------------------------------------------------------------
# Arrays to store the data that will be plotted
# FV3 data
ps    = np.zeros((N,N,6,len(modes)),dtype=np.float32)
ps_av = np.zeros((N,N,6,len(modes)),dtype=np.float32)
te    = np.zeros((N,N,6,nplots,len(modes)),dtype=np.float32)
Te = np.zeros((nplots,len(modes)))

# This loop over tiles computes the maximum errors
for k, filepath in enumerate(filepaths):
    for tile in range(0,6):
        print(filepath, tile)
        # Files to be opened
        atmos_file = filepath+"atmos_daily.tile"+str(tile+1)+".nc"
        grid_file  = filepath+"grid_spec.tile"+str(tile+1)+".nc"

        # Load the files
        data = xr.open_dataset(atmos_file, decode_times=False)
        grid = xr.open_dataset(grid_file , decode_times=False)
        area = grid['area'].values

        # Variable to be plotted
        t = nplots-1
        ps   [:,:,tile,k] = data['ps'][t,:,:].values
        ps_av[:,:,tile,k] = data['ps_av'][t,:,:].values

        for t in range(0, nplots):
           te[:,:,tile,t,k] = data['te'][t,:,:].values
           Te[t, k] =  np.sum(te[:,:,tile,t,k]*area)
           #print(t,Te[t, k], np.amin(area), np.amax(area))

#exit()
# fields to be plotted
fields = [ps,ps_av]
fnames = ['ps', 'ps_av']
titles = ['Surface pressure (Pa)', 'Time averaged surface pressure (Pa)']

# fields to be plotted
#fields = [ps_av,]
#fnames = ['ps_av',]
#fdifs = [ps_av_dif]

# plot
##############################################################################################
for field, fname, title in zip(fields, fnames, titles):
    field_32 = field[:,:,:,0]
    field_64 = field[:,:,:,1]

    fmin = min(np.amin(field_32), np.amin(field_64))
    fmax = max(np.amax(field_32), np.amax(field_64))

    # difference
    dif = abs(field_64-field_32)#/field_64
    difabs = np.amax(dif)
    auxtitle = testname+'_diff_'+fname+'_t'+str(Tf)+'_'+basename
    filename = outdir+auxtitle
    auxtitle = 'ABS(SINGLE-DOUBLE): '+title+ ', T = '+str(Tf)+' days'
    plot_scalarfield(dif, auxtitle, filename, filepaths[1], cmap_precp, 0, difabs)

    # Compute spatial mean absolute error
    spatial_mae  = np.sum(dif)/(6*N*N)

    # Compute spatial root-mean absolute error
    spatial_rmse = np.sqrt(np.sum(dif*dif)/(6*N*N))

    if(fname=='ps_av'):
       print('MAE :', spatial_rmse)
       print('RMSE:', spatial_mae)


    ##############################################################################################