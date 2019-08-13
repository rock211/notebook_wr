import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
import pyathena as pa
from multiprocessing import Pool
from matplotlib.ticker import MultipleLocator

from mpl_toolkits import axes_grid1

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

unit = pa.set_units(muH=1.4271)
print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
volpercell = 7168.*1024*1024/(128*128*896)
vpc = volpercell

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.





simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'MHD_8pc_new' ,
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
S = (':','--','-')

k = 1
# overplot Starformation rate of three different simulations
#plt.figure(figsize=(14,8))
fig =plt.figure(figsize=(8.5,12))
for j in (3,4): #range(1,7) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    Mom_1 = []
    Mom_2 = []
    Mom_3 = []
    Mom_c = []
    Mom_u = []
    Mom_w = []
    Mom_h = []
    Mom_i = []
    Mom_cuw = []
    Mom_hi = []
    Mom_t = []
    Mom_icm = []
    #plt.figure(figsize=(8, 5)) # shrink 12,5

    for tidx in range(250, stop):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        #print(ds.derived_field_list)

        # full domain information
        #print ds.domain

        # information of grid #0
        #print ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields

        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        d = ds.read_all_data('density')*unit['density'].value # density
        #den = ds.read_all_data('density')
        #P = ds.read_all_data('pressure')
        #mass = d # density times volume per cell = mass per cell
        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        #vel = ds.read_all_data('velocity')
        #T1 = ds.read_all_data('T1')
        #vel_z = vel[:,:,:,2]
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func

        surfd = np.sum(d,axis=0)*8/64.
        ax =plt.imshow(surfd,origin='lower',norm=LogNorm(),cmap='bone_r')
        plt.clim(0.1/64.,100/64.)
        #cbar = fig.colorbar()
        plt.colorbar()
        plt.show()
