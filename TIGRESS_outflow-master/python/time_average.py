import yt
import numpy as np
import glob,os
import xarray as xr

tmin=200
tmax=500

pid='MHD_4pc_new'
mvz_dir='/tigress/changgoo/{}/MVZ/'.format(pid)
mpdf_all={}

yt.funcs.mylog.setLevel(50) 
for phase in ['warm','int','hot']:
    files = glob.glob('{}{}/{}_mvz.*.h5'.format(mvz_dir,phase,pid))
    files.sort()
    print(len(files))
    time=[]
    mpdf=[]
    for f in files:
        ds = yt.load(f)
        tMyr=ds.current_time.to('Myr').v
        if (tMyr >= tmin) & (tMyr <= tmax):
            print(tMyr)
            time.append(tMyr)
            mpdf.append(ds.data['cell_mass'].to('Msun').v)
        if tMyr > tmax: 
            break
    vbin=ds.data['x_bins'].v
    zbin=ds.data['y_bins'].v
    vc=0.5*(vbin[1:]+vbin[:-1])
    zc=0.5*(zbin[1:]+zbin[:-1])

    mpdf_xr=xr.DataArray(np.array(mpdf),coords=[time,vc,zc],dims=['time','vz','z'])

    mpdf_xr.to_netcdf('../data/mpdf_{}.nc'.format(phase))
