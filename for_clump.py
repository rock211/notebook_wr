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
from astropy.io import fits

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

simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
#labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14')
labelll = ('Cold','Unstable','Warm','Ionized','Hot')

C2 = ('darkblue','orange','goldenrod','red','firebrick')

S = (':','--','-')
hh = [0.0062,0.0062,0.0062,0.02,0.152,0.152,0.25]
ylim = [0.0062,0.0062,0.0062,0.02,0.152,0.152,0.25]
ylim_c = [0.15,0.15,0.15,0.35,10,10,10]
inflow = [0,0.005002/4,0.005002,0.005002,0.005002*2,0.005002*2,0.005002*4] # icm inflow rate
k = 0


icm_c = '#03d803ff'
jj = (4,2)
fig, axs = plt.subplots(len(jj),1, figsize=(10,12), sharex = True)

for j in jj: #range(1,7) :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    if j == 5 or j==6:
        stop = 473
    else:
        stop = 499

    Mom_1 = []
    Mom_2 = []
    Mom_3 = []
    Mom_c = []
    Mom_u = []
    Mom_w = []
    Mom_i = []
    Mom_h = []
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

        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        d = ds.read_all_data('density')*unit['density'].value*8*8*8 # density in solarmass per pc3
        hdu = fits.PrimaryHDU(d)
        hdu.writeto('F:/yeongu/clump/%s_%s.fits' % (labell[j],tidx))
        print tidx