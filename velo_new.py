import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa

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
vpc = 7168.*1024*1024/(128*128*896) # volume per cell

simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc')  # 'MHD_8pc_new' ,
labell = ('No ICM', 'Weak', 'Strong', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
S = (':','--','-')

# overplot Starformation rate of three different simulations

for j in (1, 2):
    if j!=2:
        basedir = 'F:/yeongu/'
    else:
        basedir = 'D:/yeongu/'
    simid = simid_t[j]
    V0 = []
    V1 = []
    V5 = []
    for tidx in range(250, 499):  # time step 251, 331, 411, 501
        #plt.figure(figsize=(2.5,7))
        #surf = ('{}{}/surf/{}.{:04d}.surf.p'.format(basedir,simid,simid,tidx))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        vel = ds.read_all_data('velocity')
        vel_z = vel[:,:,:,2]

        #print rp.shape
        vz = np.mean(np.mean(vel_z,axis=1),axis=1)
        v0 = vz[0]
        v_n5 = vz[386] # -0.5 kpc
        v_n1= vz[323] # -1 kpc height
        V0.append(v0)
        V1.append(v_n1)
        V5.append(v_n5)
        print tidx
    t = range(250,499)
    plt.plot(t,V0,label='-3.5kpc')
    plt.plot(t,V1,label='-1kpc',alpha=0.8)
    plt.plot(t,V5,label='-0.5kpc',alpha=0.5)
    plt.legend()
    plt.show()