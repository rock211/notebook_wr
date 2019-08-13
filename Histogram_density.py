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

# In[9]:

# You can retrive the unit system with pa.set_units function.
# To make unit conversion easier, I use astropy's unit and constant.
# You may need to install astropy to use it
# please visit http://www.astropy.org/

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

#stop = 500

simid_t = ('RPS_8pc_ICM00','RPS_8pc_ICM0') # 'R8_8pc_metal','RPS_8pc_ICM0','RPS_8pc_ICM3'
labell = (r'No ICM','ICM0','ICM1','ICM2','ICM3') # 'NonICM',r'No ICM','ICM0',,'ICM3'
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

# overplot Starformation rate of three different simulations

for j in range(len(simid_t)) :
    basedir = 'F:/yeongu/'
    simid = simid_t[0]
    Mom_up = []
    data = []
    stop = 500
    #if j != len(simid_t):
    #    stop = 500
    #else:
    #    stop = 472
    #print stop,j
    for tidx in range(250,stop):  # time step 251, 331, 411, 501

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


        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        s_cut = 0.7
        if j !=0 :
            scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
            d[scalar > s_cut] = 0
        #vel = ds.read_all_data('velocity')
        T1 = ds.read_all_data('T1')
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func

        s_d = np.sum(d*8,axis=0) # surface density : vertical sum(volume density x height(8pc)) = surface density
        s_d = np.log10(s_d[s_d >0])
        #print s_d.min()
        mean = round(np.mean(s_d),3)
        medi = round(np.median(s_d),3)
        #std = np.std(s_d)
        #dist = lognorm([std],loc=mean)
        #print tidx
        #tt.append((s_d.min,s_d.max()))
        #bins = 199
        #print xx
        #shape, loc, scale = lognorm.fit(s_d, floc=0)
        #n, bins, patches = plt.hist(s_d,bins=xx,facecolor='blueviolet',alpha=0.7)#,log=True)
        #xbin = np.arange(-2,2,np.log10(1.1))
        xbin = np.linspace(-5,3,801)
        hist, binedges = np.histogram(s_d,bins = xbin)
        #print hist
        #print binedges
        #print binedges[0:-1]
        #print (binedges[1::]-binedges[0:-1])/2
        #fit = lognorm.pdf(xx-1,shape,loc=loc,scale=scale)
        #print bins
        #area_hist = .0
        #for ii in range(n.size):
        #    area_hist += (bins[ii + 1] - bins[ii]) *n[ii]

        #plt.plot(xx-1,fit*area_hist,c='k')
        #print bins
        x_c = (binedges[1::]+binedges[0:-1])/2
        binwidth = binedges[1]-binedges[0]

        data.append(hist) # for data making
        '''
        #print len(data)
        #plt.yscale('log')
        #print binedges[0:-1]+(binedges[1::]-binedges[0:-1])/2
        #print (binedges[1::]-binedges[0:-1])
        plt.bar(x_c,hist,align='center',edgecolor='k',linewidth=0.2,width =binwidth,alpha=0.7)
        #plt.bar(binedges[0:-1], hist,align='center',alpha=0.5,edgecolor='k',width = 0.1)
        plt.xlabel(r'$log$ SurfaceDensity$[M_{\odot}/pc^2]$')
        plt.ylabel('Number')
        plt.xlim(-2,3)
        #plt.ylim(0.5,1e4)
        plt.title('Hist_%s_%s' % (labell[j],tidx))
        ylim = 1500
        plt.ylim(0,ylim)
        #plt.xlim(0,10)
        #plt.text(7,5*1e5,r'scalar cut = 0.5')
        #plt.text(np.amin(bins)*1.05,ylim*0.9,r'Mean = %s' % mean)
        #plt.text(np.amin(bins) *1.05, ylim * 0.85, r'Median = %s' % medi)
        #plt.text(np.amin(bins) *1.05, ylim * 0.8, r's_cut = 0.75')
        #plt.text(np.amax(bins) * 0.75, ylim* 0.8, r'bin width = 1')
        #plt.savefig('D:/yeongu/plots/surfD_hist/Hist_%s_%s.png' % (labell[j], tidx),dpi=100)

        #plt.show()
        #plt.close()
        #print tidx
        print tidx
        '''
        print tidx
    ####### make histogram data set ###########   
    print map(str, x_c)
    row = map(str,range(250,stop))
    zdata = pd.DataFrame(data,columns=map(str, x_c),index=row)
    print zdata
    zdata.to_pickle('F:/yeongu/%s/surf_hist/surf_distribution_%s_3' % (simid,simid))
    #print zdata
    #print tt.min(),tt.max()
