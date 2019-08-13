import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from scipy.stats import lognorm

from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl


# In[8]:

# ## Unit system
#
# The unit system we choose for this simulation is
# * [length] = pc
# * [velocity] = km/s
# * [density] = 1.4271*m_h/cm^3

unit = pa.set_units(muH=1.4271)
print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
vpc = 7168.*1024*1024/(128*128*896) # volume per cell


# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

simid_t = ('RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'R8_8pc_metal'
labell = ('ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # 'NonICM',r'No ICM',
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

# overplot Starformation rate of three different simulations

for j in range(0,2) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j] # simid_t[j]
    if j == 4 :
        stop = 472
    elif j == 5 :
        stop = 401
    else :
        stop = 500
    ICM_v = []
    ICM_m = []
    for tidx in range(250,stop):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        #print(ds.derived_field_list)

        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
        icm_m_frac = scalar*d/np.sum(d) # mass fraction
        icm_v_frac = scalar/(896*128*128.) # volume fraction
        #d_icm = d[scalar > 0.7]
        '''
        ############# ICM histo ###########
        #hist, binedges = np.histogram(icm, bins=np.arange(0.7, 1, 0.01),weights=d_icm)
        hist, binedges = np.histogram(icm, bins=np.arange(0.7, 1, 0.01))
        x_c = (binedges[1::]+binedges[0:-1])/2
        binwidth = binedges[1]-binedges[0]

        plt.bar(x_c, (hist), align='center', edgecolor='k', linewidth=0.2, width=binwidth, alpha=0.7,log=1)
        plt.xlabel(r'ISM               $\longleftrightarrow$               ICM')

        plt.ylabel('Number')  # no weighted
        plt.title(r'Passive Scalar_%s %s' % (labell[j],tidx)) # No weighted
        plt.ylim(1e3,1e8) # No weighted

        #plt.ylabel('')  # weighted
        #plt.title(r'Weighted_Passive Scalar_%s %s' % (labell[j], tidx)) # Weighted
        #plt.ylim(1e-2,1e2) # weighted

        plt.show()
        #plt.savefig('D:/yeongu/plots/ICM/ICM_hist_%s_%s.png' % (labell[j], tidx),dpi=150)
        plt.close()
        print tidx
        ###################################
        '''
        print tidx
        ICM_v.append(np.sum(icm_v_frac))
        ICM_m.append(np.sum(icm_m_frac))
    xx = range(250,stop)
    plt.plot(xx,ICM_v,label='Volume Fraction')
    plt.plot(xx,ICM_m,label='Mass Fraction')
    plt.xlabel('Shot N')
    plt.ylabel('Fraction')
    plt.title(r'Passive Scalar_%s' % labell[j])
    plt.ylim(0,1)
    #plt.ylabel('') # Weighted
    #plt.title(r'Weighted_Passive Scalar_%s' % labell[j]) # Weighted
    plt.legend(loc=0)
    plt.savefig('D:/yeongu/plots/ICM/ICM_hist_%s.png' % (labell[j]), dpi=100)
    #plt.show()
    plt.close()

    print labell[j]